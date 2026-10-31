#pragma once

#include "bnp/pricing/pricing.hpp"
#include "common.cpp"
#include <Highs.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <lp_data/HConst.h>
#include <numeric>
#include <queue>
#include <vector>

int pricing_bnb_nodes = 0;
int pricing_bnb_runs = 0;
struct epbnb_state {
	bitvec nodes;	   // nodes taken
	bitvec untakeable; // nodes discarded or taken or with uncompatible colors
	double score;	   // score achieved until now
	epbnb_state(const pricing_problem &p) : nodes(p.g.n()), untakeable(p.g.n()), score(0) {}
	void take(const int u, const pricing_problem &p) {
		score -= p.cost[u];
		for(auto &[v, w] : p.g.g[u])
			if(nodes[v]) score += w;
		nodes.set1(u);
		untakeable.set1(u);
		for(const auto v : ~untakeable)
			if(!compatible(p.g.c[u], p.g.c[v])) untakeable.set1(v);
	}
	void forbid(const int u) { untakeable.set1(u); }
	bool can_branch(const pricing_problem &p) const {
		for(const auto u : nodes)
			for(const auto &[v, w] : p.g.g[u])
				if(!untakeable[v]) return 1;
		return 0;
	}
	int branch_node(const pricing_problem &p) const {
		int minv = p.g.n();
		double bv = -std::numeric_limits<double>::infinity();
		for(auto u : ~untakeable) {
			double val = -p.cost[u];
			bool adj = 0;
			std::vector<double> bc(p.g.C(), 0);
			for(const auto &[v, w] : p.g.g[u]) {
				if(nodes[v]) {
					adj = 1;
					val += w;
				} else if(!untakeable[v]) {
					bc[p.g.c[v][0]] = std::max(bc[p.g.c[v][0]], w-p.cost[v]);
				}
			}
			for(double &x : bc) val += x;
			if(adj && val>bv) {
				bv = val;
				minv = u;
			}
		}
		assert(minv!=p.g.n());
		return minv;
	}
	std::array<epbnb_state, 2> branch(const pricing_problem &p, const int minv) const {
		epbnb_state take_node(*this);
		epbnb_state forbid_node(*this);
		take_node.take(minv, p);
		forbid_node.forbid(minv);
		return {take_node, forbid_node};
	}
	void update_untakeable_reach(const pricing_problem &p) {
		if(untakeable.count_ones()==0) return;
		bitvec reach(p.g.n());
		std::vector<int> bc(p.g.C(), 0);
		std::queue<std::array<int, 2>> q;
		int cl = p.g.C(); // colors left to take
		/*{
			bitvec cvis(p.g.C());
			cvis.zero();
			for(int i = 0; i<p.g.n(); ++i)
				if(!untakeable[i] && !cvis[p.g.c[i][0]]) {
					cl += 1;
					cvis.set1(p.g.c[i][0]);
				}
		}*/
		for(auto i : nodes) cl -= p.g.c[i].size();
		for(auto u : nodes) {
			for(auto &[v, w] : p.g.g[u])
				if(!untakeable[v] && !reach[v] && w>0) { // w<=0 means it's useless in pricing preprocessing
					q.push({v, 1});
					reach.set1(v);
				}
		}
		while(!q.empty()) {
			auto [u, d] = q.front();
			q.pop();
			if(d==cl) continue;
			for(auto &[v, w] : p.g.g[u])
				if(!untakeable[v] && !reach[v] && w>0) {
					q.push({v, d+1});
					reach.set1(v);
				}
		}
		for(auto u : reach) {
			fill(bc.begin(), bc.end(), 0);
			int wsum = 0;
			for(auto &[v, w] : p.g.g[u]) {
				if(nodes[v] || (reach[v] && w>bc[p.g.c[v][0]])) {
					wsum -= bc[p.g.c[v][0]];
					bc[p.g.c[v][0]] = w;
					wsum += bc[p.g.c[v][0]];
					if(wsum>=p.cost[u]) break;
				}
			}
			if(wsum<p.cost[u]) reach.set0(u);
		}
		untakeable |= ~reach;
	}
	double upper_bound(const pricing_problem &p, const std::vector<int> &msigma) {
		// return upper_bound_surrogate_lp(p);
		update_untakeable_reach(p);
		std::vector<double> bbc(p.g.C(), 0); // best reachable value by color
		std::vector<double> bc(p.g.C(), 0);	 // best connection by color, for internal loop
		for(auto u : ~untakeable) {
			fill(bc.begin(), bc.end(), 0);
			for(auto &[v, w] : p.g.g[u]) {
				if(nodes[v]) {
					bc[p.g.c[v][0]] = w;
				} else if(!untakeable[v]) {
					double valu = msigma[u]+p.cost[u]; // -sigma[u]
					double valv = msigma[v]+p.cost[v]; // -sigma[v]
					const double muv = std::min(valu, valv);
					valu -= muv;
					valv -= muv;
					const double alpha = (1.0+valu)/(2.0+valu+valv);
					bc[p.g.c[v][0]] = std::max(bc[p.g.c[v][0]], w*alpha);
				}
			}
			double wsum = 0;
			for(const auto &sc : bc) wsum += sc;
			const double us = wsum-p.cost[u];
			if(us>bbc[p.g.c[u][0]]) bbc[p.g.c[u][0]] = us;
		}
		double ans = score;
		for(const auto &sc : bbc) ans += sc;

		return ans;
	}
	// linearize the thing
	double upper_bound_lp(const pricing_problem &p) {
		const double inf = kHighsInf;
		update_untakeable_reach(p);
		Highs h;
		h.setOptionValue("output_flag", false);
		h.changeObjectiveSense(ObjSense::kMaximize);
		std::vector<int> tk;
		std::vector<int> rtk(p.g.n(), -1);
		for(int u = 0; u<p.g.n(); ++u)
			if(!untakeable[u]) {
				rtk[u] = tk.size();
				tk.push_back(u);
			}
		if(tk.size()==0) return score;
		int coli = 0;
		for(int i = 0; i<(int)tk.size(); ++i) {
			h.addVar(0, 1);
			++coli;
		}
		for(int k = 0; k<p.g.C(); ++k) {
			vector<int> ind;
			vector<double> val;
			for(int i = 0; i<int(tk.size()); ++i) {
				int u = tk[i];
				bool contains = 0;
				for(auto &ki : p.g.c[u])
					if(ki==k) {
						contains = 1;
						break;
					}
				if(contains) {
					ind.push_back(i);
					val.push_back(1);
				}
			}
			if(ind.size()) h.addRow(-inf, 1, ind.size(), ind.data(), val.data());
		}
		for(int i = 0; i<(int)tk.size(); ++i) {
			int u = tk[i];
			for(auto [v, w] : p.g.g[u]) {
				int j = rtk[v];
				if(j>i && w>0) { // j>i => j!=-1
					h.addVar(0, 1);
					// h.changeColIntegrality(coli, HighsVarType::kImplicitInteger);
					int nz[2] = {i, coli};
					double val[2] = {-1, 1};
					h.addRow(-inf, 0, 2, nz, val);
					nz[0] = j;
					h.addRow(-inf, 0, 2, nz, val);
					h.changeColCost(coli, w);
					++coli;
				}
			}
		}
		for(int i = 0; i<(int)tk.size(); ++i) {
			int u = tk[i];
			double val = -p.cost[u];
			for(auto [v, w] : p.g.g[u]) {
				if(nodes[v]) val += w;
			}
			h.changeColCost(i, val);
		}
		h.run();
		assert(h.getModelStatus()==HighsModelStatus::kOptimal);
		return h.getObjectiveValue()+score;
	}
	// solve the surrogate with lp
	double upper_bound_surrogate_lp(const pricing_problem &p) {
		const double inf = kHighsInf;
		update_untakeable_reach(p);
		std::vector<int> tk;
		std::vector<int> rtk(p.g.n(), -1);
		std::vector<std::array<int, 3>> te;
		vector<bool> cvis(p.g.C());
		vector<int> rtc(p.g.C());
		int CC = 0;
		for(int u = 0; u<p.g.n(); ++u)
			if(!untakeable[u]) {
				rtk[u] = tk.size();
				tk.push_back(u);
				if(!cvis[p.g.c[u][0]]) {
					rtc[p.g.c[u][0]] = CC;
					CC += 1;
					cvis[p.g.c[u][0]] = 1;
				}
				for(auto &[v, w] : p.g.g[u])
					if(v>u && !untakeable[v] && w>0) { te.push_back({u, v, w}); }
			}
		if(tk.size()==0) return score;
		// alpha [0,1]   : |E'|      = te.size()
		// beta  [0,inf]    : |V' x C'| = tk.size()*CC
		// gamma [0,inf] : |C'|      = CC
		// beta[u][c[v]]>=w(u,v)*alpha[uv] and beta[v][c[u]]>=w(u,v)*(1-alpha[uv])
		// forall u:c[u]=k : gamma[k]>=val[u]+sum(for k notin c[u])(beta[u][k])
		// min sum gamma

		Highs h;
		h.setOptionValue("output_flag", false);
		h.setOptionValue("log_to_console", false);
		h.setOptionValue("solver", "simplex");
		h.setOptionValue("presolve", "off");
		h.changeObjectiveSense(ObjSense::kMinimize);
		for(int i = 0; i<(int)te.size(); ++i) { // alpha
			h.addVar(0, 1);
		}
		for(int i = 0; i<(int)tk.size()*CC; ++i) { // beta
			h.addVar(0, inf);
		}
		for(int i = 0; i<CC; ++i) { // gamma
			h.addVar(0, inf);
		}
		for(int ai = 0; ai<int(te.size()); ++ai) {
			const auto &[u, v, w] = te[ai];
			{
				// beta[u][c[v]]>=w(u,v)*alpha[uv]
				// +w(u,v)*alpha[uv]-beta[u][c[v]]<=0
				int i = rtk[u];
				int ind[2] = {ai, int(te.size()+CC*i+rtc[p.g.c[v][0]])};
				double vals[2] = {double(w), -1.0};
				h.addRow(-inf, 0, 2, ind, vals);
			}

			{
				// beta[v][c[u]]>=w(u,v)*(1-alpha[uv])
				// -alpha[uv] -beta[v][c[u]]<=-w(u,v)
				int j = rtk[v];
				int ind[2] = {ai, int(te.size()+CC*j+rtc[p.g.c[u][0]])};
				double vals[2] = {-1.0, -1.0};
				h.addRow(-inf, -w, 2, ind, vals);
			}
		}
		for(int i = 0; i<(int)tk.size(); ++i) {
			int u = tk[i];
			double valu = -p.cost[u];
			for(auto [v, w] : p.g.g[u]) {
				if(nodes[v]) valu += w;
			}
			// c[u]=k : gamma[k]>=val[u]+sum(for k notin c[u])(beta[u][k])
			// +sum(for k notin c[u])(beta[u][k]) -gamma[k]<=-val[u]
			vector<int> ind;
			vector<double> vals;
			ind.reserve(CC+1);
			vals.reserve(CC+1);
			for(int ki = 0; ki<CC; ++ki) {
				ind.push_back(te.size()+CC*i+ki);
				vals.push_back(1.0);
			}
			ind.push_back(te.size()+tk.size()*CC+rtc[p.g.c[u][0]]);
			vals.push_back(-1.0);
			h.addRow(-inf, -valu, ind.size(), ind.data(), vals.data());
		}
		for(int gi = int(te.size()+tk.size()*CC); gi<int(te.size()+tk.size()*CC+CC); ++gi) { // gammas
			h.changeColCost(gi, 1.0);
		}
		h.run();
		assert(h.getModelStatus()==HighsModelStatus::kOptimal);
		return h.getObjectiveValue()+score;
	}
	void update_primal_greedy(const pricing_problem &p, double &primal, bitvec &ans) const {
		bitvec nnodes = nodes;
		bitvec nuntakeable = untakeable;
		double nscore = score;
		bool updt;
		do {
			updt = 0;
			for(int u = 0; u<p.g.n(); ++u)
				if(!nuntakeable[u]) {
					double valu = -p.cost[u];
					for(auto [v, w] : p.g.g[u])
						if(nnodes[v]) valu += w;
					if(valu>0) {
						nnodes.set1(u);
						nuntakeable.set1(u);
						for(int v = 0; v<p.g.n(); ++v)
							if(!nuntakeable[v] && !compatible(p.g.c[u], p.g.c[v])) nuntakeable.set1(v);
						nscore += valu;
						updt = 1;
					}
				}
		} while(updt);
		for(int u = 0; u<p.g.n(); ++u)
			if(nnodes[u]) {
				double valu = -p.cost[u];
				for(auto [v, w] : p.g.g[u])
					if(nnodes[v]) valu += w;
				if(valu<0) {
					nnodes.set0(u);
					nscore -= valu;
				}
			}
		if(nscore>primal) {
			primal = nscore;
			ans = nnodes;
		}
	}
	void debug_print(const pricing_problem &p) {
		std::cerr<<"epbnb_state debug:"<<std::endl;
		std::cerr<<"score="<<score<<std::endl;
		for(int u = 0; u<p.g.n(); ++u)
			if(!untakeable[u]) {
				double valu = -p.cost[u];
				for(auto [v, w] : p.g.g[u]) {
					if(nodes[v]) valu += w;
				}
				std::cerr<<"cost("<<u<<")="<<p.cost[u]<<'\t';
				std::cerr<<"val("<<u<<")="<<valu<<"\t";
				std::cerr<<"col("<<u<<")=";
				for(auto &k : p.g.c[u]) std::cerr<<k<<',';
				std::cerr<<std::endl;
			}
		for(int u = 0; u<p.g.n(); ++u)
			if(!untakeable[u]) {
				for(const auto &[v, w] : p.g.g[u])
					if(v>u && !untakeable[v]) { std::cerr<<u<<'-'<<v<<": "<<w<<std::endl; }
			}
	}
};

void epbnb(const pricing_problem &p, const epbnb_state state, std::vector<int> &msigma, const double dual,
		   double &primal, bitvec &ans) {
	pricing_bnb_nodes += 1;
	if(state.score>primal) {
		primal = state.score;
		ans = state.nodes;
	}
	assert(state.score<=dual+0.0001);
	if(primal>=dual || !state.can_branch(p)) {
		// state.update_primal_greedy(p, primal, ans);
		return;
	}
	int u = state.branch_node(p);
	auto [s1, s2] = state.branch(p, u);
	for(const auto &[v, w] : p.g.g[u]) msigma[v] -= w;
	const double d1 = std::min(dual, s1.upper_bound(p, msigma)); // take
	for(const auto &[v, w] : p.g.g[u]) msigma[v] += w;
	const double d2 = std::min(dual, s2.upper_bound(p, msigma)); // forbid
	if(d2>d1) {
		epbnb(p, s2, msigma, d2, primal, ans);
		for(const auto &[v, w] : p.g.g[u]) msigma[v] -= w;
		epbnb(p, s1, msigma, d1, primal, ans);
		for(const auto &[v, w] : p.g.g[u]) msigma[v] += w;
	} else {
		for(const auto &[v, w] : p.g.g[u]) msigma[v] -= w;
		epbnb(p, s1, msigma, d1, primal, ans);
		for(const auto &[v, w] : p.g.g[u]) msigma[v] += w;
		epbnb(p, s2, msigma, d2, primal, ans);
	}
}
void epbnb0(const pricing_problem &p, epbnb_state state, double &primal, bitvec &ans) {
	pricing_bnb_nodes += 1;
	std::vector<int> ord(p.g.n(), 0);
	std::iota(ord.begin(), ord.end(), 0);
	std::vector<int> msigma(p.g.n(), 0);
	// std::sort(ord.begin(),ord.end(),[&](int i, int j){return p.cost[i]<p.cost[j];});
	for(auto &u : ord) {
		epbnb_state s1 = state;
		s1.take(u, p);
		for(const auto &[v, w] : p.g.g[u]) msigma[v] -= w;
		const double d1 = s1.upper_bound(p, msigma);
		epbnb(p, s1, msigma, d1, primal, ans);
		state.forbid(u);
		for(const auto &[v, w] : p.g.g[u]) msigma[v] += w;
	}
}

std::vector<mcolumn> exact_pricing_bnb(pricing_problem p) {
	pricing_bnb_runs += 1;
	epbnb_state s0(p);
	double primal = 0;
	bitvec ans = s0.nodes;
	epbnb0(p, s0, primal, ans);
	if(primal>1e-6) {
		mcolumn c;
		double val = primal;
		for(int i = 0; i<p.g.n(); ++i)
			if(ans[i]) {
				c.nodes.push_back(i);
				val += p.cost[i];
			}
		c.value = int(round(val));
		return {c};
	}
	return {};
}
