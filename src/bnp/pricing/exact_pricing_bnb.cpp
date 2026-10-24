#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
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
	/*
		split a node with multiple colors evenly across those colors
		upper_bound2 kept updated across epbnb_states
		this is hopelessly bugged
	*/
	double upper_bound2(const pricing_problem& p) {
		struct epbnb_dual2 {
			const double eps = 0.01;
			struct edge {
				int u;
				int v;
				int w;
				double alpha; // 1 is towards u, 0 is towards v
				double val(int i) const {
					if(i==u) return alpha*u;
					else return (1-alpha)*v;
				}
			};

			// color associated to each node
			std::vector<int> col;
			// list of edges
			std::vector<edge> edges;
			// for each node, for each color, connections in decreasing order of value
			std::vector<std::vector<std::vector<int>>> g;
			// for each node, total value
			std::vector<double> value;
			// base node value, without the edges
			std::vector<double> base_value;
			// for each color, nodes sorted by decreasing value
			std::vector<std::vector<int>> bbc;

			int n() const {
				return int(col.size());
			}
			int C() const {
				return g[0].size();
			}

			int can_reduce_edge(int ei, int u) const {
				int v = edges[ei].u;
				if(v==u) v = edges[ei].v;
				if(edges[ei].alpha+eps>=1.0 && edges[ei].u==v) return 0;
				if(edges[ei].alpha-eps<=0.0 && edges[ei].v==v) return 0;
				if(edges[g[v][col[u]][0]].val(v) > eps + edges[ei].val(v)) return 1; // can increase bc edge is not best for v for u's color
				if(value[bbc[col[v]][0]] > eps + value[v]) return 2; // can increase bc v is not best of its color
				return 0;
			}
			std::vector<std::array<int,3>> can_reduce_node(int u) const {
				for(auto &gc: g[u]) if(!gc.empty()) {
					std::vector<std::array<int,3>> to_reduce;
					int eii=0;
					for(eii=0; eii<(int)gc.size(); ++eii) {
						if(edges[gc[eii]].val(u)+eps < edges[gc[0]].val(u)) return to_reduce;
						int t = can_reduce_edge(gc[eii], u);
						if(t==0) break;
						to_reduce.push_back({gc[eii],u,t-1});
					}
					if(eii==(int)gc.size()) return to_reduce;
				}
				return {};
			}
			std::vector<std::array<int,3>> can_reduce_color(int c) const {
				std::vector<std::array<int,3>> to_reduce;
				for(int ui=0; ui<(int)bbc[c].size(); ++ui) {
					if(value[ui]+eps<value[bbc[c][0]]) return to_reduce;
					auto t = can_reduce_node(bbc[c][ui]);
					if(t.empty()) return {};
					to_reduce.insert(to_reduce.end(), std::make_move_iterator(t.begin()), std::make_move_iterator(t.end()));
				}
				return to_reduce;
			}
			void reduce_edge(int ei, int u, int div) {
				int v = edges[ei].u;
				if(v==u) v = edges[ei].v;
				double diff = 0.0;
				if(edges[g[v][col[u]][0]].val(v) > eps + edges[ei].val(v)) {
					diff = edges[g[v][col[u]][0]].val(v) - edges[ei].val(v);
				} else if(value[bbc[col[v]][0]] > eps + value[v]) {
					assert(div!=0);
					diff = (value[bbc[col[v]][0]] - value[v])/(double)div;
				}
				diff /= edges[ei].w;
				if(edges[ei].u==v) edges[ei].alpha += diff;
				if(edges[ei].alpha>1.0) edges[ei].alpha = 1.0;
				return;
			}
			void resort() {
				for(int u=0; u<n(); ++u) {
					value[u] = base_value[u];
					for(auto&cg: g[u]) if(!cg.empty()) {
						sort(cg.begin(),cg.end(),[&](int ei, int ej){
							return edges[ei].val(u) > edges[ej].val(u);
						});
						value[u]+=edges[cg[0]].val(u);
					}
				}
				for(auto&bc: bbc) {
					sort(bc.begin(),bc.end(),[&](int u, int v){ return value[u]>value[v]; });
				}
			}
			bool try_reduce_color(int c) {
				std::vector<std::array<int,3>> to_reduce = can_reduce_color(c);
				if(to_reduce.empty()) return 0;
				std::vector<int> div(n());
				for(auto &[ei,u,t]: to_reduce) {
					int v = edges[ei].u == u ? edges[ei].v : edges[ei].u;
					div[v] += t;
				}
				for(auto &[ei,u,t]: to_reduce) {
					int v = edges[ei].u == u ? edges[ei].v : edges[ei].u;
					reduce_edge(ei, u, div[v]);
				}
				resort();
				return 1;
			}
			epbnb_dual2(const pricing_problem& p, const epbnb_state& state) {
				std::vector<int> colmap(p.g.n(),-1); // original color -> compressed color
				int coln=0;
				std::vector<int> nodemap(p.g.n(),-1); // original node -> compressed node
				int noden=0;
				for(int u=0; u<p.g.n(); ++u) if(!state.untakeable[u]) {
					nodemap[u] = noden++;
					if(colmap[p.g.c[u][0]]==-1) {
						colmap[p.g.c[u][0]]=coln++;
					}
				}
				base_value.resize(noden,0);
				value.resize(noden,0);
				g.resize(noden,std::vector<std::vector<int>>(coln));
				col.resize(noden);
				bbc.resize(coln);
				for(int u=0; u<p.g.n(); ++u) if(!state.untakeable[u]) {
					col[nodemap[u]] = colmap[p.g.c[u][0]];
					base_value[nodemap[u]] = -p.cost[u];
					for(auto &[v,w]: p.g.g[u]) if(state.nodes[v]) {
						base_value[nodemap[u]]+=w;
					}
				}
				for(int u=0; u<p.g.n(); ++u) if(!state.untakeable[u]) {
					for(auto &[v,w]: p.g.g[u]) if(v>u && !state.untakeable[v]) {
						int uu = nodemap[u];
						int vv = nodemap[v];
						g[uu][col[vv]].push_back(edges.size());
						g[vv][col[uu]].push_back(edges.size());
						edges.push_back(edge{uu,vv,w,0.5});
					}
				}
				for(int u=0; u<n(); ++u) {
					bbc[col[u]].push_back(u);
				}
				resort();
			}
			double ub() const {
				double ans = 0;
				for(int c=0; c<C(); ++c) {
					int u = bbc[c][0];
					ans += value[u];
				}
				return ans;
			}
		};
		update_untakeable_reach(p);
		epbnb_dual2 d2(p, *this);
		bool updated=0;
		std::cerr<<"AA"<<std::endl;
		do {
			updated=0;
			for(int i=0; i<d2.C(); ++i) {
				bool up = d2.try_reduce_color(i);
				if(up) std::cerr<<"Updated "<<d2.ub()<<std::endl;
				if(up) updated=1;
			}
		} while(updated);
		std::cerr<<"BB"<<std::endl;
		return d2.ub();
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
    const double d1 = std::min(dual,s1.upper_bound2(p));
    const double d2 = std::min(dual,s2.upper_bound2(p));
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
		const double d1 = s1.upper_bound2(p);
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

