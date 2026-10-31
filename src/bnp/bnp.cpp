#pragma once

#include "bnp/master.hpp"
#include "common.cpp"
#include "heuristic/heuristic_karger.cpp"
#include "mgraph.cpp"
#include "scorer.cpp"
#include "bnp/hotstart.cpp"
#include <iostream>
#include <map>
#include <queue>
#include <set>

int bnp_nodes = 0;
/// global set of columns, only insert sorted vectors of nodes with indeces in the original graph
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

struct BnpNode {
	mgraph mg;
	int base; // base value (sum of edges already merged in mg)
	int dual; // dual value computed by the parent node
	int depth;
	bool operator<(const BnpNode &o) const { return std::make_pair(dual, -depth)<std::make_pair(o.dual, -o.depth); }
};
void bnp_bfs(const mgraph mg0, int &primal, std::vector<std::vector<int>> &best) {
	if(mg0.is_trivial()) { // base case: all components are colorful
		int pv = mg0.m();
		if(pv>primal) {
			primal = pv;
			auto comps = mg0.split_components();
			best = std::vector<std::vector<int>>();
			for(auto &comp : comps) { best.push_back(comp.orig_nodes()); }
		}
		return;
	}
	BnpNode n0 = BnpNode{mg0, 0, mg0.m(), 0};
	std::priority_queue<BnpNode> maxheap;
	maxheap.push(n0);
	while(!maxheap.empty()) {
		bnp_nodes += 1;
		auto [mg, base, dual, depth] = maxheap.top();
		maxheap.pop();
		std::cerr<<"BNPNODE "<<bnp_nodes<<": depth="<<depth<<" base="<<base<<" dual="<<dual
				 <<" primal="<<primal<<std::endl;
		if(primal>=dual) {
			std::cerr<<"prune1"<<std::endl;
			continue;
		}

		Master m(mg, columns_for(mg));
		m.setup_highs();
		double dv = base+m.run(primal-base);
		for(auto col : m.get_columns()) { columns.insert(col); }
		const std::vector<double> cv = m.highs.getSolution().col_value;
		assert(cv.size()==m.columns.size());
		dual = std::min(dual, int(dv+1e-4));
		if(primal>=dual) {
			std::cerr<<"prune2"<<std::endl;
			continue;
		}

		auto [pv, pcols] = m.primal_ilp(primal);
		pv += base;
		if(pv>primal) {
			primal = pv;
			best = pcols;
			std::cerr<<"NEW PRIMAL: "<<primal<<std::endl;
		}
		if(primal>=dual) {
			std::cerr<<"prune3"<<std::endl;
			continue;
		}

		// find where to split
		std::array<int, 2> split = {-1, -1};
		{
			std::map<std::array<int, 2>, double> tmp;
			for(int i = 0; i<(int)cv.size(); ++i)
				if(cv[i]>0.0001) {
					sort(m.columns[i].nodes.begin(), m.columns[i].nodes.end()); // to call intersect
					for(auto &u : m.columns[i].nodes) {
						for(auto &v : intersect(mg.g[u], m.columns[i].nodes, u)) { tmp[{u, v}] += cv[i]; }
					}
				}
			double md = 2;
			for(auto &[uv, w] : tmp) {
				if(fabs(w-0.5)<md) {
					md = fabs(w-0.5);
					split = uv;
				}
			}
		}
		assert(split[0]!=-1);
		std::cerr<<"branch on ("<<split[0]<<","<<split[1]<<")"<<std::endl;
		// TODO: split components
		{ // merge edge
			mgraph mg1(mg);
			int w = mg1.merge_edge(split[0], split[1]);
			assert(w);
			maxheap.push({mg1, base+w, dual, depth+1});
		}

		{ // remove edge
			mgraph mg2(mg);
			assert(mg2.remove_edge(split[0], split[1]));
			maxheap.push({mg2, base, dual, depth+1});
		}
	}
}

/// mg already is a single component
void bnp(const mgraph &mg, const int depth, const int base, int dual, int &primal,
		 std::vector<std::vector<int>> &best) {
	bnp_nodes += 1;
	if(mg.is_trivial()) { // base case: all components are colorful
		int pv = base+mg.m();
		if(pv>primal) {
			primal = pv;
			auto comps = mg.split_components();
			best = std::vector<std::vector<int>>();
			for(auto &comp : comps) { best.push_back(comp.orig_nodes()); }
		}
		return;
	}
	std::cerr<<"BNPNODE "<<bnp_nodes<<": depth="<<depth<<" base="<<base<<" dual="<<dual
			 <<" primal="<<primal<<std::endl;

	Master m(mg, columns_for(mg));
	m.setup_highs();
	double dv = base+m.run(primal-base);
	for(auto col : m.get_columns()) { columns.insert(col); }
	const std::vector<double> cv = m.highs.getSolution().col_value;
	assert(cv.size()==m.columns.size());
	dual = std::min(dual, int(dv+1e-4));
	if(primal>=dual) {
		std::cerr<<"prune1"<<std::endl;
		return;
	}
	auto [pv, pcols] = m.primal_ilp(primal-base);
	pv += base;
	if(pv>primal) {
		primal = pv;
		best = pcols;
		std::cerr<<"NEW PRIMAL: "<<primal<<std::endl;
	}
	if(primal>=dual) {
		std::cerr<<"prune2"<<std::endl;
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
	bnp(mg1, depth+1, base+w, dual, primal, best);

	if(primal>=dual) {
		std::cerr<<"prune second child\n";
		return;
	}

	mgraph mg2(mg);
	assert(mg2.remove_edge(split[0], split[1]));
	bnp(mg2, depth+1, base, dual, primal, best);
}

std::vector<std::vector<int>> bnp_main(problem p, bool hot_start = 1) {
	columns = std::set<std::vector<int>>(); // reset columns
	// single-element columns for feasibility are not needed
	// hot start columns
	if(hot_start) {
		for(auto col : gen_hotstart_columns(p)) columns.insert(col);
	}
	mgraph mg0(p);

	std::vector<std::vector<int>> ans;
	auto ccs = mg0.split_components();
	int sol = 0;
	for(auto &cc : ccs) {
		std::cerr<<"component: n="<<cc.n()<<" partitions_bound="<<cc.tight_partition_bound()<<std::endl;
		std::vector<std::vector<int>> ansi;
		int primal = -1;
		bnp(cc, 0, 0, cc.m(), primal, ansi);
		sol += primal;
		std::move(ansi.begin(), ansi.end(), std::back_inserter(ans));
	}
	assert(score(p, ans)==sol);
	return ans;
}
