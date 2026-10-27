#pragma once

#include "bnp/pricing/pricing.hpp"
#include "common.cpp"
#include <Highs.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <lp_data/HConst.h>
#include <numeric>
#include <queue>
#include <vector>

// TODO: optimize
int pricing_bnb_nodes = 0;
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
		for(int v = 0; v<p.g.n(); ++v)
			if(!untakeable[v]) {
				if(!compatible(p.g.c[u], p.g.c[v])) untakeable.set1(v);
			}
	}
	void forbid(const int u) { untakeable.set1(u); }
	std::array<epbnb_state, 2> branch(const pricing_problem &p) const {
		int minv = p.g.n();
		for(int u = 0; u<p.g.n() && minv==p.g.n(); ++u)
			if(nodes[u]) {
				for(auto &[v, w] : p.g.g[u])
					if(!untakeable[v]) {
						minv = v;
						break;
					}
			}

		epbnb_state take_node(*this);
		epbnb_state forbid_node(*this);
		take_node.take(minv, p);
		forbid_node.forbid(minv);
		return {take_node, forbid_node};
	}
	void update_untakeable_reach(const pricing_problem &p) {
		bitvec reach(p.g.n());
		bitvec cvis(p.g.C());
		std::vector<int> bc(p.g.C(), 0);
		bool updated = 0;
		do {
			int cl = 0; // colors left to take
			{
				cvis.zero();
				for(int i = 0; i<p.g.n(); ++i)
					if(!untakeable[i] && !cvis[p.g.c[i][0]]) {
						cl += 1;
						cvis.set1(p.g.c[i][0]);
					}
			}
			updated = 0;
			reach.zero();
			std::queue<std::array<int, 2>> q;
			for(int u = 0; u<p.g.n(); ++u)
				if(nodes[u]) {
					for(auto &[v, w] : p.g.g[u])
						if(!untakeable[v] && !reach[v]) {
							q.push({v, 1});
							reach.set1(v);
						}
				}
			while(!q.empty()) {
				auto [u, d] = q.front();
				q.pop();
				if(d==cl) continue;
				for(auto &[v, w] : p.g.g[u])
					if(!untakeable[v] && !reach[v]) {
						q.push({v, d+1});
						reach.set1(v);
					}
			}
			for(int u = 0; u<p.g.n(); ++u)
				if(reach[u]) {
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
			for(int u = 0; u<p.g.n(); ++u)
				if(!untakeable[u] && !reach[u]) {
					updated = 1;
					untakeable.set1(u);
				}
		} while(updated);
	}
	double upper_bound(const pricing_problem &p) {
		update_untakeable_reach(p);
		std::vector<double> bbc(p.g.C(), 0); // best reachable value by color
		std::vector<double> bc(p.g.C(), 0);	 // best connection by color, for internal loop
		// TODO: better bound
		// TODO: for each color save best even and best odd so that you can floor after dividing sum of wsum by 2(?)
		for(int u = 0; u<p.g.n(); ++u)
			if(!untakeable[u]) {
				fill(bc.begin(), bc.end(), 0);
				for(auto &[v, w] : p.g.g[u]) {
					if(nodes[v]) {
						bc[p.g.c[v][0]] = w;
					} else if(!untakeable[v]) {
						const double val = (4+p.cost[v])/(8.0+p.cost[u]+p.cost[v]);
						// const double val = 0.5;
						bc[p.g.c[v][0]] = std::max(bc[p.g.c[v][0]], w*val /*/2.0*/);
					}
				}
				double wsum = 0;
				for(auto &sc : bc) wsum += sc;
				const double us = wsum-p.cost[u];
				// std::cerr<<"aoeu "<<u<<' '<<wsum<<' '<<us<<std::endl;
				if(us>bbc[p.g.c[u][0]]) bbc[p.g.c[u][0]] = us;
			}
		double ans = score;
		for(auto &sc : bbc) ans += sc;

		// const double dlp = upper_bound_lp(p);
		const double dslp = upper_bound_surrogate_lp(p);
		// assert(dslp<=dlp);
		if(dslp>ans+0.0001) {
			std::cerr<<"dslp="<<dslp<<" ans="<<ans<<std::endl;
			debug_print(p);
			upper_bound_surrogate_lp(p,1);
		}
		assert(dslp<=ans+0.0001);
		return dslp;
		// return ans;
	}
	// linearize the thing
	double upper_bound_lp(const pricing_problem &p) {
		const double inf = kHighsInf;
		update_untakeable_reach(p);
		Highs h;
		// h.setOptionValue("seed",42); // TODO seed
		h.setOptionValue("output_flag", false);
		h.changeObjectiveSense(ObjSense::kMaximize);
		std::vector<int> tk;
		std::vector<int> rtk(p.g.n(), -1);
		for(int u = 0; u<p.g.n(); ++u) if(!untakeable[u]) {
			rtk[u]=tk.size();
			tk.push_back(u);
		}

		int coli = 0;
		for(int i=0; i<(int)tk.size(); ++i) {
			h.addVar(0, 1);
			++coli;
		}
		for(int i = 0; i<int(tk.size()); ++i) {
			for(int j = i+1; j<int(tk.size()); ++j) {
				int u = tk[i];
				int v = tk[j];
				if(!compatible(p.g.c[u], p.g.c[v])) {
					int ind[2] = {i, j};
					double val[2] = {1, 1};
					h.addRow(-inf, 1, 2, ind, val);
				}
			}
		}
		if(tk.size()==0) return score;
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
		/*
		if(h.getModelStatus()!=HighsModelStatus::kOptimal) {
			debug_print(p);
			std::cerr<<"lp "<<h.modelStatusToString(h.getModelStatus())<<std::endl;
		}
		*/
		assert(h.getModelStatus()==HighsModelStatus::kOptimal);
		return h.getObjectiveValue()+score;
	}
	// solve the surrogate with lp
	double upper_bound_surrogate_lp(const pricing_problem &p, bool verbose=0) {
		const double inf = kHighsInf;
		update_untakeable_reach(p);
		std::vector<int> tk;
		std::vector<int> rtk(p.g.n(), -1);
		std::vector<std::array<int,3>> te;
		vector<bool> cvis(p.g.C());
		vector<int> rtc(p.g.C());
		int CC = 0;
		for(int u = 0; u<p.g.n(); ++u) if(!untakeable[u]) {
			rtk[u]=tk.size();
			tk.push_back(u);
			if(!cvis[p.g.c[u][0]]) {
				rtc[p.g.c[u][0]]=CC;
				CC+=1;
				cvis[p.g.c[u][0]]=1;
			}
			for(auto &[v, w]: p.g.g[u]) if(v>u && !untakeable[v] && w>0) {
				te.push_back({u,v,w});
			}
		}
		if(tk.size()==0) return score;
		// alpha [0,1]   : |E'|      = te.size()
		// beta  [0,inf]    : |V' x C'| = tk.size() * CC
		// gamma [0,inf] : |C'|      = CC
		// beta[u][c[v]] >= w(u,v)*alpha[uv] and beta[v][c[u]] >= w(u,v)*(1-alpha[uv])
		// forall u:c[u]=k : gamma[k] >= val[u] + sum(for k notin c[u])(beta[u][k])
		// min sum gamma

		Highs h;
		// h.setOptionValue("seed",42); // TODO seed
		h.setOptionValue("output_flag", false);
		h.changeObjectiveSense(ObjSense::kMinimize);
		for(int i=0; i<(int)te.size(); ++i) { // alpha
			h.addVar(0,1);
		}
		for(int i=0; i<(int)tk.size()*CC; ++i) { // beta
			h.addVar(0, inf);
		}
		for(int i=0; i<CC; ++i) { // gamma
			h.addVar(0, inf);
		}
		for(int ai=0; ai<int(te.size()); ++ai) {
			const auto &[u,v,w] = te[ai];
			{
				// beta[u][c[v]] >= w(u,v)*alpha[uv]
				// +w(u,v)*alpha[uv] - beta[u][c[v]] <= 0
				if(verbose) {
					std::cerr<<"beta["<<u<<"]["<<rtc[p.g.c[v][0]]<<"] >= "<<w<<" * alpha["<<u<<","<<v<<"]\n";
				}
				int i = rtk[u];
				int ind[2] = {ai, int(te.size() + CC*i + rtc[p.g.c[v][0]])};
				double vals[2] = {double(w),-1.0};
				h.addRow(-inf,0,2,ind,vals);
			}
			
			{
				// beta[v][c[u]] >= w(u,v)*(1-alpha[uv])
				// -alpha[uv] -beta[v][c[u]] <= -w(u,v)
				if(verbose) {
					std::cerr<<"beta["<<v<<"]["<<rtc[p.g.c[u][0]]<<"] >= "<<w<<" * (1-alpha["<<u<<","<<v<<"])\n";
				}
				int j = rtk[v];
				int ind[2] = {ai, int(te.size() + CC*j + rtc[p.g.c[u][0]])};
				double vals[2] = {-1.0,-1.0};
				h.addRow(-inf, -w, 2, ind, vals);
			}
		}
		for(int i=0; i<(int)tk.size(); ++i) {
			int u = tk[i];
			double valu = -p.cost[u];
			for(auto [v, w] : p.g.g[u]) {
				if(nodes[v]) valu += w;
			}
			// c[u]=k : gamma[k] >= val[u] + sum(for k notin c[u])(beta[u][k])
			// +sum(for k notin c[u])(beta[u][k]) -gamma[k] <= -val[u]
			vector<int> ind;
			vector<double> vals;
			ind.reserve(CC+1);
			vals.reserve(CC+1);
			if(verbose) std::cerr<<"gamma["<<rtc[p.g.c[u][0]]<<"] >= "<<valu;
			for(int ki=0; ki<CC; ++ki) {
				ind.push_back(te.size() + CC*i + ki);
				vals.push_back(1.0);
				if(verbose) std::cerr<<" +beta["<<u<<"]["<<ki<<"]";
			}
			if(verbose) std::cerr<<std::endl;
			ind.push_back(te.size() + tk.size()*CC + rtc[p.g.c[u][0]]);
			vals.push_back(-1.0);
			h.addRow(-inf, -valu, ind.size(), ind.data(), vals.data());
		}
		for(int gi=int(te.size()+tk.size()*CC); gi<int(te.size()+tk.size()*CC+CC); ++gi) { // gammas
			h.changeColCost(gi, 1.0);
		}
		h.run();
		/*
		if(h.getModelStatus()!=HighsModelStatus::kOptimal) {
			debug_print(p);
			std::cerr<<"slp "<<h.modelStatusToString(h.getModelStatus())<<std::endl;
		}
		*/
		assert(h.getModelStatus()==HighsModelStatus::kOptimal);
		if(verbose) {
			const auto& sol = h.getSolution().col_value;
			for(int i=0; i<(int)te.size(); ++i) { // alpha
				std::cerr<<"alpha: ("<<te[i][0]<<','<<te[i][1]<<','<<te[i][2]<<"):"<<sol[i]<<std::endl;
			}
			for(int i=0; i<(int)tk.size()*CC; ++i) { // beta
				std::cerr<<"beta["<<tk[i/CC]<<"]["<<i%CC<<"]="<<sol[i+te.size()]<<std::endl;
			}
			for(int i=0; i<CC; ++i) { // gamma
				std::cerr<<"gamma["<<i<<"]="<<sol[te.size()+tk.size()*CC+i]<<std::endl;
			}
			std::cerr<<"objvalue="<<h.getObjectiveValue()<<std::endl;
		}
		return h.getObjectiveValue()+score;
	}
	void debug_print(const pricing_problem& p) {
		std::cerr<<"epbnb_state debug:"<<std::endl;
		std::cerr<<"score="<<score<<std::endl;
		for(int u = 0; u<p.g.n(); ++u) if(!untakeable[u]) {
				double valu = -p.cost[u];
				for(auto [v, w] : p.g.g[u]) {
					if(nodes[v]) valu += w;
				}
				std::cerr<<"cost("<<u<<")="<<p.cost[u]<<'\t';
				std::cerr<<"val("<<u<<")="<<valu<<"\t";
				std::cerr<<"col("<<u<<")=";
				for(auto &k:p.g.c[u]) std::cerr<<k<<',';
				std::cerr<<std::endl;
		}
		for(int u = 0; u<p.g.n(); ++u) if(!untakeable[u]) {
			for(const auto &[v,w]: p.g.g[u]) if(v>u && !untakeable[v]) {
				std::cerr<<u<<'-'<<v<<": "<<w<<std::endl;
			}
		}
	}
};

void epbnb(const pricing_problem &p, const epbnb_state state, const double dual, double &primal, bitvec &ans) {
	pricing_bnb_nodes += 1;
	if(state.score>primal) {
		primal = state.score;
		ans = state.nodes;
	}
	if(state.score>dual+0.0001) return;
	assert(state.score<=dual+0.0001);
	if(primal>=dual) return;
	auto [s1, s2] = state.branch(p);
	const double d1 = std::min(dual, s1.upper_bound(p));
	const double d2 = std::min(dual, s2.upper_bound(p));
	if(true || d2>d1 /*s2.score>s1.score*/) {
		if(d2>primal) { // TODO: resimplify
			epbnb(p, s2, d2, primal, ans);
			if(primal > d2+0.0001) {
				s2.debug_print(p);
				std::cerr<<"primal="<<primal<<std::endl;
				std::cerr<<"dual2="<<d2<<std::endl;
				s2.upper_bound_surrogate_lp(p,1);
				exit(1);
			}
		}
		if(primal > dual+0.0001) return;
		if(d1>primal) {
			epbnb(p, s1, d1, primal, ans);
			if(primal > d1+0.0001) {
				s1.debug_print(p);
				std::cerr<<"primal="<<primal<<std::endl;
				std::cerr<<"dual1="<<d1<<std::endl;
				s1.upper_bound_surrogate_lp(p,1);
				exit(1);
			}
		}
		if(primal > dual+0.0001) return;
	} else {
		epbnb(p, s1, d1, primal, ans);
		epbnb(p, s2, d2, primal, ans);
	}
}
void epbnb0(const pricing_problem &p, epbnb_state state, double &primal, bitvec &ans) {
	pricing_bnb_nodes += 1;
	std::vector<int> ord(p.g.n(), 0);
	std::iota(ord.begin(), ord.end(), 0);
	// std::sort(ord.begin(),ord.end(),[&](int i, int j){return p.cost[i]<p.cost[j];});
	for(auto &i : ord) {
		epbnb_state s1 = state;
		s1.take(i, p);
		const double d1 = s1.upper_bound(p);
		epbnb(p, s1, d1, primal, ans);
		state.forbid(i);
	}
}

std::vector<mcolumn> exact_pricing_bnb(pricing_problem p) {
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
