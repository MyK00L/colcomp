#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>
#include <queue>
#include <array>
#include "common.cpp"
#include "bnp/pricing/pricing.hpp"

// TODO: optimize
int nnodes=0;
struct epbnb_state {
	bitvec nodes; // nodes taken
	bitvec untakeable; // nodes discarded or taken or with uncompatible colors
	double score; // score achieved until now
	epbnb_state(const pricing_problem& p): nodes(p.g.n()), untakeable(p.g.n()), score(0) {}
	void take(const int u, const pricing_problem& p) {
		score-=p.cost[u];
		for(auto &[v,w]: p.g.g[u]) if(nodes[v]) score+=w;
		nodes.set1(u);
		untakeable.set1(u);
		for(int v=0; v<p.g.n(); ++v) if(!untakeable[v]) {
			if(!compatible(p.g.c[u],p.g.c[v])) untakeable.set1(v);
		}
	}
	void forbid(const int u) {
		untakeable.set1(u);
	}
	std::array<epbnb_state,2> branch(const pricing_problem& p) const {
		int minv=p.g.n(); 
		for(int u=0; u<p.g.n()&&minv==p.g.n(); ++u) if(nodes[u]) {
			for(auto &[v,w]:p.g.g[u]) if(!untakeable[v]) {
				minv = v;
				break;
			}
		}
		/*
		   vector<int> ord(p.g.n(),0);
		   std::iota(ord.begin(),ord.end(),0);
		   std::sort(ord.begin(),ord.end(),[&](int i, int j){return p.cost[i]<p.cost[j];});
		   vector<bool> reach(p.g.n());
		   for(auto &u:ord) if(nodes[u]) {
		   for(auto &[v,w]:p.g.g[u]) if(!untakeable[v]) {
		   reach[v]=1;
		   }
		   }
		   for(auto &u:ord) if(reach[u]) {
		   minv=u;
		   break;
		   }
		   */
		epbnb_state take_node(*this);
		epbnb_state forbid_node(*this);
		take_node.take(minv, p);
		forbid_node.forbid(minv);
		return {take_node,forbid_node};
	}
	void update_untakeable_reach(const pricing_problem& p) {
		bitvec reach(p.g.n());
		bitvec cvis(p.g.C());
		std::vector<int> bc(p.g.C(), 0);
		bool updated=0;
		do {
			int cl = 0; // colors left to take
			{
				cvis.zero();
				for(int i=0; i<p.g.n(); ++i) if(!untakeable[i] && !cvis[p.g.c[i][0]]) {
					cl+=1;
					cvis.set1(p.g.c[i][0]);
				}
			}
			updated=0;
			reach.zero();
			std::queue<std::array<int,2> > q;
			for(int u=0; u<p.g.n(); ++u) if(nodes[u]) {
				for(auto &[v,w]:p.g.g[u]) if(!untakeable[v] && !reach[v]) {
					q.push({v,1});
					reach.set1(v);
				}
			}
			while(!q.empty()) {
				auto [u,d] = q.front();
				q.pop();
				if(d==cl) continue;
				for(auto &[v,w]:p.g.g[u]) if(!untakeable[v] && !reach[v]) {
					q.push({v,d+1});
					reach.set1(v);
				}
			}
			for(int u=0; u<p.g.n(); ++u) if(reach[u]) {
				fill(bc.begin(),bc.end(),0);
				int wsum=0;
				for(auto &[v,w]: p.g.g[u]) {
					if(nodes[v] || (reach[v] && w > bc[p.g.c[v][0]] ) ) {
						wsum -= bc[p.g.c[v][0]];
						bc[p.g.c[v][0]] = w;
						wsum += bc[p.g.c[v][0]];
						if(wsum>=p.cost[u]) break;
					}
				}
				if(wsum<p.cost[u]) reach.set0(u);
			}
			for(int u=0; u<p.g.n(); ++u) if(!untakeable[u] && !reach[u]) {
				updated=1;
				untakeable.set1(u);
			}
		} while(updated);
	}
	double upper_bound(const pricing_problem& p) {
		++nnodes;
		update_untakeable_reach(p);
		std::vector<double> bbc(p.g.C(),0); // best reachable value by color
        std::vector<double> bc(p.g.C(), 0); // best connection by color, for internal loop
		// TODO: better bound
		// TODO: for each color save best even and best odd so that you can floor after dividing sum of wsum by 2(?)
		for(int u=0; u<p.g.n(); ++u) if(!untakeable[u]) {
           	fill(bc.begin(),bc.end(),0);
           	for(auto &[v,w]: p.g.g[u]) {
               	if(nodes[v]) {
                   	bc[p.g.c[v][0]] = w;
               	} else if (!untakeable[v]) {
                   	const double val = (4+p.cost[v])/(8.0+p.cost[u]+p.cost[v]);
					// const double val = 0.5;
					bc[p.g.c[v][0]] = std::max(bc[p.g.c[v][0]],w*val/*/2.0*/);
				}
           	}
			double wsum=0;
           	for(auto &sc:bc) wsum+=sc;
			const double us=wsum - p.cost[u];
			if(us>bbc[p.g.c[u][0]]) bbc[p.g.c[u][0]]=us;
		}
        double ans = score;
        for(auto &sc: bbc) ans+=sc;
		return ans;
    }
};

void epbnb(const pricing_problem& p, const epbnb_state state, const double dual, double& primal, bitvec& ans) {
    if(state.score>primal) {
        primal=state.score;
        ans=state.nodes;
    }
	assert(state.score<=dual+0.0001);
    if(primal>=dual) return;
	auto [s1,s2] = state.branch(p);
    const double d1 = std::min(dual,s1.upper_bound(p));
    const double d2 = std::min(dual,s2.upper_bound(p));
    if(d2>d1/*s2.score>s1.score*/) {
        epbnb(p,s2,d2,primal,ans);
        epbnb(p,s1,d1,primal,ans);
    } else {
        epbnb(p,s1,d1,primal,ans);
		epbnb(p,s2,d2,primal,ans);
    }
}
void epbnb0(const pricing_problem& p, epbnb_state state, const double dual, double& primal, bitvec& ans) {
	std::vector<int> ord(p.g.n(),0);
	std::iota(ord.begin(),ord.end(),0);
	// std::sort(ord.begin(),ord.end(),[&](int i, int j){return p.cost[i]<p.cost[j];});
	for(auto &i:ord) {
		epbnb_state s1 = state;
		s1.take(i,p);
		const double d1 = s1.upper_bound(p);
		epbnb(p,s1,d1,primal,ans);
		state.forbid(i);
	}
}

std::vector<mcolumn> exact_pricing_bnb(pricing_problem p) {
	// nnodes=0;
	epbnb_state s0(p);
	double primal=0;
    double dual = p.g.orig_n*p.g.orig_n;
    bitvec ans = s0.nodes;
    epbnb0(p,s0,dual,primal,ans);
    // std::cerr<<"nnodes="<<nnodes<<"\n";
    if(primal>1e-6) {
        mcolumn c;
        double val = primal;
        for(int i=0; i<p.g.n(); ++i) if(ans[i]) {
            c.nodes.push_back(i);
            val+=p.cost[i];
        }
        c.value = int(round(val));
		return {c};
    }
    return {};
}

std::vector<mcolumn> hot_start_cols(mgraph mg) {
    // TODO
    return {};
}

