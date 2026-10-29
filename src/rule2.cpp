#include "common.cpp"
#include <algorithm>
#include <queue>

// TODO: generate candidates like in 2012
// generate candidates like karger would

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
		n=0;
		int mindeg=n;
		for(auto u:gb) mp[u]=n++;
		for(auto u: gb) {
			for(const auto v: og[u]) if(gb[v] && mp[v]>mp[u]) {
				g[mp[u]].push_back({mp[v],1,int(g[mp[v]].size())});
				g[mp[v]].push_back({mp[u],1,int(g[mp[u]].size()-1)});
			}
			if(int(g[mp[u]].size())<mindeg) {
				mindeg=g[mp[u]].size();
				v0=u;
				if(mindeg<t) return 0;
			}
		}
	}
	std::vector<int> ff(n,-1);
	std::vector<Edge*> fe(n); // I hate writing flows
	bitvec S(n);
	bitvec U(n); U.flipall();
	U.set0(v0);
	S.set1(v0);
	for(const auto e: g[v0]) U.set0(e.v);
	while(U.begin()!=U.end()) {
		int v = *U.begin();
		int flow=0;
		// flow from v to S
		// flow through T
		for(auto& e1: g[v]) if(!U[e1.v]) {
			for(auto &e2: g[e1.v]) if(S[e2.v]) {
				assert(e1.w==1 && e2.w==1);
				e1.w=0;
				e2.w=0;
				flow+=1;
				// don't need to set backwards flow here
			}
		}
		// flow through U with bfs
		bool found;
		do {
			found=0;
			std::queue<int> q;
			std::fill(ff.begin(),ff.end(),-1);
			q.push(v);
			ff[v]=v;
			while(!q.empty()) {
				int u = q.front(); q.pop();
				if(S[u]) { // found a path
					while(ff[u]!=u) {
						fe[u]->w-=1;
						g[fe[u]->v][fe[u]->rwj].w+=1;
						u=ff[u];
					}
					flow+=1;
					found=1;
					break;
				}
				for(auto& e: g[u]) if(e.w>0 && ff[e.v] != -1) {
					ff[e.v] = u;
					fe[e.v] = &e;
					q.push(e.v);
				}
			}
		} while(found);

		if(flow<t) return 0;
		U.set0(v);
		S.set1(v);
		for(const auto& ei:g[v]) U.set0(ei.v);
		// reset flow
		for(auto u: ~S) for(auto& ei: g[u]) ei.w=1;
	}
	return 1;
}

int rule2(problem& p) {

}
