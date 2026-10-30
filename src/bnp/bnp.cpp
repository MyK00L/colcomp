#pragma once

#include "bnp/master.hpp"
#include "common.cpp"
#include "heuristic/heuristic_karger.cpp"
#include "mgraph.cpp"
#include "scorer.cpp"
#include <iostream>
#include <map>
#include <set>

/// global set of columns
std::set<std::vector<int>> columns;
std::vector<mcolumn> columns_for(const mgraph &mg) {
	std::vector<mcolumn> ans;
	// single-element columns for feasibility
	// for(int i=0; i<mg.n(); ++i) ans.push_back({0,{i}});

	std::vector<int> group_of_node(mg.orig_n, -1); // what group each node is in in the subproblem, -1 if finalized
	for(int i = 0; i<(int)mg.groups.size(); ++i)
		for(auto &u : mg.groups[i]) group_of_node[u] = i;
	std::vector<int> ng(mg.groups.size(), 0);
	for(const auto &col : columns) {
		bool ok = 1; // if column is compatible with subproblem
		for(auto &u : col) {
			if(group_of_node[u]<0) {
				ok = 0;
				break;
			}
			ng[group_of_node[u]] += 1;
		}
		if(ok)
			for(auto &u : col) {
				if(ng[group_of_node[u]]!=(int)mg.groups[group_of_node[u]].size()) {
					ok = 0;
					break;
				}
			}
		if(ok) { // add column
			int val = 0;
			std::vector<int> mgcol;
			for(auto &uu : col)
				if(mg.groups[group_of_node[uu]][0]==uu) {
					int u = group_of_node[uu];
					mgcol.push_back(u);
					for(auto &[v, w] : mg.g[u])
						if(v>u && ng[v]) { val += w; }
				}
			ans.push_back({val, mgcol});
		}
		// reset ng
		for(auto &u : col) {
			if(group_of_node[u]<0) break;
			ng[group_of_node[u]] = 0;
		}
	}
	return ans;
}

std::vector<int> intersect(const std::vector<std::array<int, 2>> gu, const std::vector<int> ci, const int exupto) {
	auto i = gu.begin();
	auto ie = gu.end();
	auto j = ci.begin();
	auto je = ci.end();
	while(i!=ie && i->at(0)<=exupto) ++i;
	while(j!=je && *j<=exupto) ++j;
	std::vector<int> ans;
	while(i!=ie && j!=je) {
		if(i->at(0)<*j) ++i;
		else if(*j<i->at(0)) ++j;
		else {
			ans.push_back(*j);
			++i;
			++j;
		}
	}
	return ans;
}

// TODO: limit pricing time to branch more in bnp if necessary
// (exp in dense graphs) ?

int bnp_nodes = 0;
/// mg already is a single component
void bnp(const mgraph &mg, const int base, int dual, int &primal, std::vector<std::vector<int>> &best) {
	bnp_nodes += 1;
	std::cerr<<"\n[[ "<<primal<<' '<<dual<<'\n';
	if(mg.is_trivial()) { // base case: all components are colorful
		int pv = base+mg.m();
		if(pv>primal) {
			primal = pv;
			auto comps = mg.split_components();
			best = std::vector<std::vector<int>>();
			for(auto &comp : comps) { best.push_back(comp.orig_nodes()); }
		}
		std::cerr<<"\n]]\n";
		return;
	}

	Master m(mg, columns_for(mg));
	m.setup_highs();
	double dv = base+m.run(primal-base);
	for(auto col : m.get_columns()) { columns.insert(col); }
	const std::vector<double> cv = m.highs.getSolution().col_value;
	assert(cv.size()==m.columns.size());
	dual = std::min(dual, int(dv+1e-4));
	if(primal>=dual) {
		std::cerr<<"\n]]\n";
		return;
	}
	std::cerr<<"PRIMAL:"<<primal<<" DUAL:"<<dual<<std::endl;
	auto [pv, pcols] = m.primal_ilp(primal);
	pv += base;
	if(pv>primal) {
		primal = pv;
		best = pcols;
		std::cerr<<"NEW PRIMAL: "<<primal<<std::endl;
	}
	if(primal>=dual) {
		std::cerr<<"\n]]\n";
		return;
	}

	// find where to split
	std::map<std::array<int, 2>, double> tmp;
	for(int i = 0; i<(int)cv.size(); ++i)
		if(cv[i]>0.0001) {
			sort(m.columns[i].nodes.begin(), m.columns[i].nodes.end()); // to call intersect
			for(auto &u : m.columns[i].nodes) {
				for(auto &v : intersect(mg.g[u], m.columns[i].nodes, u)) { tmp[{u, v}] += cv[i]; }
			}
		}
	double md = 2;
	std::array<int, 2> split = {-1, -1};
	for(auto &[uv, w] : tmp) {
		if(fabs(w-0.5)<md) {
			md = fabs(w-0.5);
			split = uv;
		}
	}
	assert(split[0]!=-1);
	std::cerr<<"branch on ("<<split[0]<<","<<split[1]<<")"<<std::endl;
	// TODO: split components
	mgraph mg1(mg);
	int w = mg1.merge_edge(split[0], split[1]);
	assert(w);
	bnp(mg1, base+w, dual, primal, best);

	if(primal>=dual) {
		std::cerr<<"\n]]\n";
		return;
	}

	mgraph mg2(mg);
	assert(mg2.remove_edge(split[0], split[1]));
	bnp(mg2, base, dual, primal, best);

	std::cerr<<"\n]]\n";
}

std::vector<std::vector<int>> bnp_main(problem p, bool hot_start = 1) {
	// single-element columns for feasibility
	columns = std::set<std::vector<int>>();
	for(int i = 0; i<(int)p.n(); ++i) columns.insert({i});
	// hot start columns
	if(hot_start) {
		for(auto col : heuristic_karger(p)) columns.insert(col);
	}
	mgraph mg0(p);

	std::vector<std::vector<int>> ans;
	auto ccs = mg0.split_components();
	int sol = 0;
	for(auto &cc : ccs) {
		std::cerr<<"component: n="<<cc.n()<<" partitions_bound="<<cc.tight_partition_bound()<<std::endl;
		std::vector<std::vector<int>> ansi;
		int primal = -1;
		bnp(cc, 0, cc.m(), primal, ansi);
		sol += primal;
		std::move(ansi.begin(), ansi.end(), std::back_inserter(ans));
	}
	assert(score(p, ans)==sol);
	return ans;
}
