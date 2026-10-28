#pragma once

#include "bnp/pricing/pricing.hpp"
#include "common.cpp"
#include <Highs.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <lp_data/HConst.h>
#include <numeric>
#include <queue>
#include <vector>

#include "bnp/pricing/exact_pricing_bnb.cpp"

// int pricing_bnb_nodes = 0;
struct epbnb_state2 {
	Highs h;
	bitvec nodes;	   // nodes taken
	bitvec untakeable; // nodes discarded or taken or with uncompatible colors
	double score;	   // score achieved until now
	struct Undo {
		int took;
		std::vector<int> untook;
	};
	std::vector<Undo> und;
	std::vector<std::array<int,3>> te;
	std::vector<double> valus;
	epbnb_state2(const pricing_problem &p) : nodes(p.g.n()), untakeable(p.g.n()), score(0), valus(p.g.n(),0.0) {
		const double inf = kHighsInf;
		for(int u = 0; u<p.g.n(); ++u) if(!untakeable[u]) {
			for(auto &[v, w]: p.g.g[u]) if(v>u && !untakeable[v] && w>0) {
				te.push_back({u,v,w});
			}
		}
		// alpha [0,1]   : |E'|      = te.size()
		// beta  [0,inf]    : |V' x C'| = p.g.n() * p.g.C()
		// gamma [0,inf] : |C'|      = p.g.C()
		// beta[u][c[v]] >= w(u,v)*alpha[uv] and beta[v][c[u]] >= w(u,v)*(1-alpha[uv])
		// forall u:c[u]=k : gamma[k] >= val[u] + sum(for k notin c[u])(beta[u][k])
		// min sum gamma
		h.clear();
		h.setOptionValue("output_flag", false);
		h.setOptionValue("log_to_console", false);
		h.setOptionValue("solver", "simplex");
		h.setOptionValue("presolve", "off");
		h.changeObjectiveSense(ObjSense::kMinimize);
		for(int i=0; i<(int)te.size(); ++i) { // alpha
			h.addVar(0,1);
		}
		for(int i=0; i<p.g.n()*p.g.C(); ++i) { // beta
			h.addVar(0, inf);
		}
		for(int i=0; i<p.g.C(); ++i) { // gamma
			h.addVar(0, inf);
		}
		for(int ai=0; ai<int(te.size()); ++ai) {
			const auto &[u,v,w] = te[ai];
			{
				// beta[u][c[v]] >= w(u,v)*alpha[uv]
				// +w(u,v)*alpha[uv] - beta[u][c[v]] <= 0
				int ind[2] = {ai, int(te.size() + p.g.C()*u + p.g.c[v][0])};
				double vals[2] = {double(w),-1.0};
				h.addRow(-inf,0,2,ind,vals);
			}
			
			{
				// beta[v][c[u]] >= w(u,v)*(1-alpha[uv])
				// -alpha[uv] -beta[v][c[u]] <= -w(u,v)
				int ind[2] = {ai, int(te.size() + p.g.C()*v + p.g.c[u][0])};
				double vals[2] = {-1.0,-1.0};
				h.addRow(-inf, -w, 2, ind, vals);
			}
		}
		for(int u=0; u<p.g.n(); ++u) {
			double valu = -p.cost[u];
			for(auto [v, w] : p.g.g[u]) {
				if(nodes[v]) valu += w;
			}
			valus[u]=valu;
			// c[u]=k : gamma[k] >= val[u] + sum(for k notin c[u])(beta[u][k])
			// +sum(for k notin c[u])(beta[u][k]) -gamma[k] <= -val[u]

			vector<int> ind;
			vector<double> vals;
			ind.reserve(p.g.C()+1);
			vals.reserve(p.g.C()+1);
			for(int ki=0; ki<p.g.C(); ++ki) {
				ind.push_back(te.size() + p.g.C()*u + ki);
				vals.push_back(1.0);
			}
			ind.push_back(te.size() + p.g.n()*p.g.C() + p.g.c[u][0]);
			vals.push_back(-1.0);
			h.addRow(-inf, -valu, ind.size(), ind.data(), vals.data());
		}
		for(int gi=int(te.size()+p.g.n()*p.g.C()); gi<int(te.size()+p.g.n()*p.g.C()+p.g.C()); ++gi) { // gammas
			h.changeColCost(gi, 1.0);
		}
	}
	void set_untakeable(int u) {
		assert(!untakeable[u]);
		const double inf = kHighsInf;
		und.back().untook.push_back(u);
		untakeable.set1(u);
		// +sum(for k notin c[u])(beta[u][k]) -gamma[k] <= -val[u]
		// vvv
		// +sum(for k notin c[u])(beta[u][k]) -gamma[k] <= inf
		int rowi = 2*te.size() + u;
		h.changeRowBounds(rowi, -inf, inf);
	}
	void unset_untakeable(int u) {
		assert(untakeable[u]);
		const double inf = kHighsInf;
		untakeable.set0(u);
		int rowi = 2*te.size() + u;
		h.changeRowBounds(rowi, -inf, -valus[u]);
	}
	void undo(const pricing_problem& p) {
		for(auto ui=und.back().untook.rbegin(); ui!=und.back().untook.rend(); ++ui) {
			const int u = *ui;
			unset_untakeable(u);
		}
		if(und.back().took!=-1) {
			const int u = und.back().took;
			for(auto &[v,w]: p.g.g[u]) if(!untakeable[v]) {
				valus[v] -= w;
				int rowi = 2*te.size() + v;
				h.changeRowBounds(rowi, -kHighsInf, -valus[v]);
			}
			score -= valus[u];
			nodes.set0(u);
		}
		und.pop_back();
	}
	void take(const int u, const pricing_problem &p) {
		assert(!untakeable[u]);
		und.push_back({u,std::vector<int>()});

		nodes.set1(u);

		score += valus[u];
		for(auto &[v,w]: p.g.g[u]) if(!untakeable[v]) {
			valus[v] += w;
			int rowi = 2*te.size() + v;
			h.changeRowBounds(rowi, -kHighsInf, -valus[v]);
		}
		set_untakeable(u);
		for(int v = 0; v<p.g.n(); ++v)
			if(!untakeable[v]) {
				if(!compatible(p.g.c[u], p.g.c[v])) {
					set_untakeable(v);
				}
			}
		update_untakeable_reach(p);
	}
	void forbid(const int u, const pricing_problem& p, bool root=0) {
		und.push_back({-1,std::vector<int>()});
		set_untakeable(u);
		if(!root) update_untakeable_reach(p);
	}
	bool can_branch(const pricing_problem& p) const {
		for(int u=0; u<p.g.n(); ++u)
			if(nodes[u])
				for(const auto&[v,w]: p.g.g[u])
					if(!untakeable[v])
						return 1;
		return 0;
	}
	int branch_node(const pricing_problem &p) const {
		int minv = p.g.n();
		double bv=-std::numeric_limits<double>::infinity();
		for(int u=0; u<p.g.n(); ++u) if(!untakeable[u]) {
			double val = -p.cost[u];
			bool adj = 0;
			std::vector<double> bc(p.g.C(),0);
			for(const auto &[v,w]: p.g.g[u]) {
				if(nodes[v]) {
					adj=1;
					val+=w;
				} else if(!untakeable[v]) {
					bc[p.g.c[v][0]] = std::max(bc[p.g.c[v][0]], w-p.cost[v]);
				}
			}
			for(double &x:bc) val+=x;
			if(adj && val>bv) {
				bv=val;
				minv = u;
			}
		}
		assert(minv!=p.g.n());
		return minv;
	}
	int branch_node_lp(const pricing_problem &p) const {
		const vector<double>& lpsol = h.getSolution().col_value;
		assert(lpsol.size() == te.size() + p.g.n()*p.g.C() + p.g.C());
		double muval = -std::numeric_limits<double>::infinity();
		int mu=-1;
		for(int u=0; u<p.g.n(); ++u) if(!untakeable[u]) {
			double uval = valus[u];
			for(int k=0; k<p.g.C(); ++k) {
				int buki = int(te.size() + p.g.C()*u + k);
				uval+=lpsol[buki];
			}
			if(uval>muval) {
				muval=uval;
				mu=u;
			}
		}
		assert(mu!=-1);
		return mu;
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
					set_untakeable(u);
				}
		} while(updated);
	}
	// solve the surrogate with lp
	double upper_bound_surrogate_lp() {
		h.run();
		assert(h.getModelStatus()==HighsModelStatus::kOptimal);
		return h.getObjectiveValue()+score;
	}
	double upper_bound_fast(const pricing_problem &p) {
		std::vector<double> bbc(p.g.C(), 0); // best reachable value by color
		std::vector<double> bc(p.g.C(), 0);	 // best connection by color, for internal loop
		
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
				if(us>bbc[p.g.c[u][0]]) bbc[p.g.c[u][0]] = us;
			}
		double ans = score;
		for(auto &sc : bbc) ans += sc;

		return ans;
	}
	double upper_bound(const pricing_problem& p, double primal) {
		double ubf = upper_bound_fast(p);
		if(ubf<=primal) return ubf;
		return upper_bound_surrogate_lp();
	}
	std::pair<double,bitvec> get_primal(const pricing_problem& p) {
		vector<std::pair<double,int> > bbc(p.g.C(),std::make_pair(-std::numeric_limits<double>::infinity(),-1));
		const vector<double>& lpsol = h.getSolution().col_value;
		assert(lpsol.size() == te.size() + p.g.n()*p.g.C() + p.g.C());
		for(int u=0; u<p.g.n(); ++u) if(!untakeable[u]) {
			double uval = valus[u];
			for(int k=0; k<p.g.C(); ++k) {
				int buki = int(te.size() + p.g.C()*u + k);
				uval+=lpsol[buki];
			}
			bbc[p.g.c[u][0]]=std::max(bbc[p.g.c[u][0]],std::make_pair(uval,u));
		}
		sort(bbc.rbegin(),bbc.rend());
		bitvec nnodes(p.g.n());
		double nscore = score;
		for(const auto &[_,u]: bbc) if(u!=-1) {
			nscore+=valus[u];
			for(auto &[v,w]: p.g.g[u]) if(nnodes[v]) nscore+=w;
			nnodes.set1(u);
		}
		std::reverse(bbc.begin(),bbc.end());
		for(const auto&[_,u]: bbc) if(u!=-1) {
			double uscore = valus[u];
			for(auto &[v,w]: p.g.g[u]) if(nnodes[v]) uscore+=w;
			if(uscore<0) {
				nscore-=uscore;
				nnodes.set0(u);
			}
		}
		return {nscore, nnodes|nodes};
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

void epbnb2(const pricing_problem &p, epbnb_state2& state, double dual, double &primal, bitvec &ans) {
	pricing_bnb_nodes += 1;
	
	{
		auto [nscore,nans] = state.get_primal(p);
		if(nscore>primal) {
			primal=nscore;
			ans=nans;
		}
	}
	if(state.score>primal) {
		primal = state.score;
		ans = state.nodes;
	}
	dual = std::min(dual,state.upper_bound(p, primal));

	assert(state.score<=dual+0.0001);
	if(primal>=dual || !state.can_branch(p)) return;

	int bn = state.branch_node_lp(p);
	
	state.take(bn, p);
	epbnb2(p,state, dual,primal,ans);
	state.undo(p);
	
	if(primal>=dual) return;

	state.forbid(bn, p);
	epbnb2(p,state,dual,primal,ans);
	state.undo(p);

}
void epbnb02(const pricing_problem &p, epbnb_state2& state, double &primal, bitvec &ans) {
	pricing_bnb_nodes += 1;
	std::vector<int> ord(p.g.n(), 0);
	std::iota(ord.begin(), ord.end(), 0);
	// std::sort(ord.begin(),ord.end(),[&](int i, int j){return p.cost[i]<p.cost[j];});
	for(auto &i : ord) {
		state.take(i, p);
		const double d1 = state.upper_bound_surrogate_lp();
		epbnb2(p, state, d1, primal, ans);
		state.undo(p);
		state.forbid(i, p, 1);
	}
}

std::vector<mcolumn> exact_pricing_bnb2(pricing_problem p) {
	epbnb_state2 s0(p);
	double primal = 0;
	bitvec ans = s0.nodes;
	epbnb02(p, s0, primal, ans);
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
