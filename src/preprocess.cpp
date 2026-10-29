#pragma once

#include "common.cpp"
#include "heuristic/heuristic_karger.cpp"
#include <algorithm>
#include <iostream>
#include <queue>

// TODO: optimize stuff

bool check_candidate(const problem& p, const std::vector<std::vector<int>>& og, const int C, const bitvec& gb) {
	int t=0;
	bitvec cols(C);
	for(auto u: gb) {
		for(const auto v: og[u]) if(!gb[v]) {
			t+=1;
			cols.set1(p.color[v]);
		}
	}
	int n = gb.count_ones();
	if(t >= n) return 0; // can't be t-edge connected
	for(auto u: gb) cols.set0(p.color[u]);
	if(cols.count_ones()) return 0; // a color in N is not present in G_B

	// check for t-edge connectivity with Matula 1987 EC algo
	struct Edge{
		int v;
		int w;
		int rwj;
	};
	std::vector<std::vector<Edge>> g(n);
	int v0=0;
	{
		std::vector<int> mp(p.n(), -1);
		int mindeg=n;
		n=0;
		for(auto u:gb) mp[u]=n++;
		for(auto u: gb) {
			for(const auto v: og[u]) if(gb[v] && mp[v]>mp[u]) {
				g[mp[u]].push_back({mp[v],1,int(g[mp[v]].size())});
				g[mp[v]].push_back({mp[u],1,int(g[mp[u]].size()-1)});
			}
			if(int(g[mp[u]].size())<mindeg) {
				mindeg=g[mp[u]].size();
				v0=mp[u];
				if(mindeg<t) return 0;
			}
		}
	}
	std::vector<int> ff(n,-1);
	std::vector<std::array<int*,2> > fe(n); // I hate writing flows
	bitvec S(n);
	bitvec U(n); U.flipall();
	U.set0(v0); // oh no v0 >= n
	S.set1(v0);
	for(const auto e: g[v0]) U.set0(e.v);
	while(U.begin()!=U.end()) {
		int v = *U.begin();
		int flow=0;
		// flow from v to S
		// flow through nt
		for(auto vt: g[v]) if(!U[vt.v]) {
			for(auto vs: g[vt.v]) if(S[vs.v]) {
				flow+=1;
				vt.w-=1;
				vs.w-=1;
				// no need to add backflow for source and sink
			}
		}
		// flow through nu
		bool found=1;
		while(found && flow<t) {
			found=0;
			std::queue<int> q;
			std::fill(ff.begin(),ff.end(),-1);
			q.push(v);
			ff[v]=v;
			while(!q.empty() && flow<t) {
				int u = q.front(); q.pop();
				if(S[u]) { // found a path
					while(u!=ff[u]) {
						assert(*fe[u][0] > 0);
						*fe[u][0]-=1;
						*fe[u][1]+=1;
						u=ff[u];
					}
					found = 1;
					flow += 1;
					break;
				}
				for(auto& e: g[u]) if(e.w>0 && ff[e.v] == -1) {
					ff[e.v] = u;
					fe[e.v] = {&e.w, &g[e.v][e.rwj].w};
					assert(*fe[e.v][0]>0);
					q.push(e.v);
				}
			}
		}

		if(flow<t) return 0;
		U.set0(v);
		S.set1(v);
		for(const auto& ei:g[v]) U.set0(ei.v);
		// reset flow
		for(auto u: ~S) for(auto& ei: g[u]) ei.w=1;
	}
	return 1;
}

bool check_candidate_mine(const problem& p, const std::vector<std::vector<int>>& og, const int C, const bitvec& gb) {
	int t=0;
	std::vector<int> colus(C, -1);
	for(auto u: gb) colus[p.color[u]] = u;
	for(auto u: gb) {
		for(const auto v: og[u]) if(!gb[v]) {
			t+=1;
			if(colus[p.color[v]]==-1) return 0; // a color in N is not present in G_B
		}
	}
	int n = gb.count_ones();
	if(t >= n) return 0; // all cuts are < t
	struct Edge{
		int v;
		int w;
		int rwj;
	};
	std::vector<std::vector<Edge>> g(n);
	std::vector<std::array<int,2>> tocheck;
	{
		std::vector<int> mp(p.n(), -1);
		n=0;
		for(auto u:gb) mp[u]=n++;
		for(auto u: gb) {
			for(const auto v: og[u]) if(gb[v] && mp[v]>mp[u]) {
				g[mp[u]].push_back({mp[v],1,int(g[mp[v]].size())});
				g[mp[v]].push_back({mp[u],1,int(g[mp[u]].size()-1)});
			} else if(!gb[v]) {
				const int o = mp[colus[p.color[v]]];
				if(mp[u]<o) tocheck.push_back({mp[u], o});
				else tocheck.push_back({o,mp[u]});
			}
		}
		sort(tocheck.begin(),tocheck.end());
		tocheck.erase( unique( tocheck.begin(), tocheck.end() ), tocheck.end() );
	}
	std::vector<int> ff(n,-1);
	std::vector<std::array<int*,2> > fe(n); // I hate writing flows
	for(const auto [u0,u1]: tocheck) {
		int flow=0;
		bool found=1;
		while(found && flow<t) {
			found=0;
			std::queue<int> q;
			std::fill(ff.begin(),ff.end(),-1);
			q.push(u0);
			ff[u0]=u0;
			while(!q.empty()) {
				int u = q.front(); q.pop();
				if(u==u1) {
					while(u!=ff[u]) {
						assert(*fe[u][0] > 0);
						*fe[u][0]-=1;
						*fe[u][1]+=1;
						u=ff[u];
					}
					found = 1;
					flow += 1;
					break;
				}
				for(auto& e: g[u]) if(e.w>0 && ff[e.v] == -1) {
					ff[e.v] = u;
					fe[e.v] = {&e.w, &g[e.v][e.rwj].w};
					q.push(e.v);
				}
			}
		}
		if(flow<t) return 0;
		// reset flow
		for(int u=0; u<(int)g.size(); ++u) for(auto& ei: g[u]) ei.w=1;
	}
	return 1;
}

/// In g: remove all edges in or adjacent to sg
/// In p: remove all edges adjacent to sg, but not inside
void remove_subgraph(problem& p, std::vector<std::vector<int>>& g, const bitvec& sg) {
	for(int u=0; u<(int)g.size(); ++u) {
		if(sg[u]) {
			g[u].clear();
		} else {
			for(int ei=0; ei<(int)g[u].size(); ++ei) {
				if(sg[g[u][ei]]) {
					std::swap(g[u][ei],g[u].back());
					g[u].pop_back();
					--ei;
				}
			}
		}
	}
	for(int ei=0; ei<(int)p.edges.size(); ++ei) {
		const auto [u,v] = p.edges[ei];
		if(sg[u]^sg[v]) {
			std::swap(p.edges[ei],p.edges.back());
			p.edges.pop_back();
			--ei;
		}
	}
}
void rule2(problem& p, std::vector<std::vector<int>> g, const int C, bool use_mine) {
	int nc1=0,nc2=0;
	for(auto comp: heuristic_karger(p)) {
		bitvec nodes(p.n());
		for(const auto u:comp) nodes.set1(u);
		bool c1 = check_candidate(p, g, C, nodes);
		bool c2 = use_mine?check_candidate_mine(p, g, C, nodes):c1;
		if(c1) assert(c2);
		if(c2) {
			nc2+=c2;
			nc1+=c1;
			remove_subgraph(p, g, nodes);
			break;
		}
	}
	for(int u0=0; u0<int(p.n()); ++u0) if(g[u0].size()) { // otherwise, it has been removed
		bitvec nodes(p.n());
		bitvec neigh(p.n());
		bitvec cols(C);
		nodes.set1(u0);
		cols.set1(p.color[u0]);
		for(int v: g[u0]) neigh.set1(v);
		for(;;) {
			int minn = p.n()*2;
			int minu = -1;
			for(const auto u : neigh) if(!cols[p.color[u]]) {
				int unn = 0;
				for(const auto v: g[u]) {
					if(nodes[v]) unn-=1;
					else unn+=1;
				}
				if(unn < minn) {
					minn = unn;
					minu = u;
				}
			}
			if(minu==-1) break;
			cols.set1(p.color[minu]);
			nodes.set1(minu);
			neigh.set0(minu);
			for(const auto v: g[minu]) if(!nodes[v]) neigh.set1(v);
			bool c1 = check_candidate(p, g, C, nodes);
			bool c2 = use_mine?check_candidate_mine(p, g, C, nodes):c1;
			if(c1) assert(c2);
			if(c2) {
				nc2+=c2;
				nc1+=c1;
				remove_subgraph(p, g, nodes);
				break;
			}
		}
	}
	std::cerr<<nc1<<' '<<nc2<<std::endl;
}

// remove colorful components
void rule1(problem& p, std::vector<std::vector<int>> g, const int C) {
	int nc=0;
	std::vector<bool> vis(p.n());
	for(int u0=0; u0<(int)p.n(); ++u0) if(!vis[u0]) {
		bool good=1;
		std::vector<bool> col(C);
		std::queue<int> q;
		bitvec comp(p.n());
		q.push(u0);
		while(!q.empty()) {
			int u = q.front(); q.pop();
			comp.set1(u);
			if(col[p.color[u]]) good=0;
			col[p.color[u]] = 1;
			for(const int v: g[u]) if(!vis[v]) {
				vis[v]=1;
				q.push(v);
			}
		}
		if(good) {
			remove_subgraph(p, g, comp);
			nc+=1;
		}
	}
	std::cerr<<"rule1: "<<nc<<std::endl;
}
// TODO: rule3?

void preprocess(problem& p, bool use_mine=0) {
	std::vector<std::vector<int>> g(p.n());
	for(auto &[u,v]: p.edges) {
		g[u].push_back(v);
		g[v].push_back(u);
	}
	int C = *std::max_element(p.color.begin(),p.color.end())+1;
	rule1(p,g,C);
	rule2(p,g,C,use_mine);
	sort(p.edges.begin(),p.edges.end());
}
