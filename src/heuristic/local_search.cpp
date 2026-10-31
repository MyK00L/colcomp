#include <algorithm>
#include <queue>
#include <vector>

#include "common.cpp"

std::vector<int> detect_cycle(const std::vector<std::vector<std::array<int,2>>>& g, const std::vector<int> pre) {
	std::vector<int> stk;
	std::vector<bool> instk(g.size());
	std::vector<bool> vis(g.size());
	for(int u0=0; u0<(int)g.size(); ++u0) if(!vis[u0]) {
		for(int u=u0; u!=-1; u=pre[u]) {
			if(!vis[u]) {
				vis[u]=1;
				stk.push_back(u);
				instk[u]=1;
			} else {
				if(instk[u]) {
					std::vector<int> cycle;
					while(stk.back()!=u) {
						cycle.push_back(stk.back());
						stk.pop_back();
					}
					cycle.push_back(u);
					return cycle;
				}
				break;
			}
		}
		for(int u: stk) instk[u]=0;
		stk.clear();
	}
	return {};
}

// https://konaeakira.github.io/posts/using-the-shortest-path-faster-algorithm-to-find-negative-cycles.html
std::vector<int> find_negative_cycle_spfa(const std::vector<std::vector<std::array<int,2>>>& g) {
	std::vector<int> len(g.size()), dist(g.size()), pre(g.size(),-1);
	std::queue<int> q;
	std::vector<bool> inq(g.size(),1);
	for(int i=0; i<(int)g.size(); ++i) q.push(i);
	int it=0;
	while(!q.empty()) {
		int u = q.front(); q.pop();
		inq[u]=0;
		for(const auto& [v,w]: g[u]) if(dist[u]+w<dist[v]) {
			it+=1;
			pre[v]=u;
			dist[v]=dist[u]+w;
			if(!inq[v]) {
				q.push(v);
				inq[v]=1;
			}
			if(it==(int)g.size()) {
				it=0;
				const auto cycle = detect_cycle(g,pre);
				if(cycle.size()) return cycle;
			}
		}
	}
	return detect_cycle(g,pre);
}

// large scale neighbourhood search with improvement graph/cyclic exchanges

int local_search_color(const problem& p, const std::vector<std::vector<int>> og, std::vector<int>& compof, const int color, std::vector<std::vector<int>>& comps) {
	std::vector<std::vector<std::array<int,2>>> g(comps.size()+1);
	std::vector<int> tocomp(comps.size());
	// edge weights are negated bc I want to find negative cycles
	for(auto &cmp: comps) for(auto &u: cmp) if((int)p.color[u]==color) {
		for(auto &v: og[u]) tocomp[compof[v]]+=1;
		for(auto &v: og[u]) if(compof[v]!=compof[u] && tocomp[compof[v]]!=0) {
			g[compof[u]].push_back({compof[v],-tocomp[compof[v]]+tocomp[compof[u]]});
			tocomp[compof[v]]=0;
		}
		g[g.size()-1].push_back({compof[u],0});
		g[compof[u]].push_back({int(g.size()-1),tocomp[compof[u]]});
		tocomp[compof[u]]=0; // reset tocomp for next cycle
		break; // only one node with a color for component
	}
	auto cycle = find_negative_cycle_spfa(g);
	int imp=0;
	for(int i=0; i<(int)cycle.size(); ++i) {
		int u = cycle[i];
		int v = cycle[(i+1)%cycle.size()];
		bool found=0;
		for(auto &[vv,w]: g[u]) if(vv==v) { // assuming there is no double edges, which should be true
			imp+=w;
			found=1;
			break;
		}
		assert(found);
	}
	assert(imp<0 || cycle.empty());
	if(cycle.size()) {
		comps.push_back({});
		int u = -1;
		for(int cii=0; cii<(int)comps[cycle.back()].size(); ++cii) if((int)p.color[comps[cycle.back()][cii]]==color) {
			std::swap(comps[cycle.back()][cii],comps[cycle.back()].back());
			u = comps[cycle.back()].back();
			comps[cycle.back()].pop_back();
			compof[u]=-1;
			break;
		}
		for(auto ci: cycle) {
			bool found = 0;
			for(int cii=0; cii<(int)comps[ci].size(); ++cii) if((int)p.color[comps[ci][cii]]==color) {
				std::swap(comps[ci][cii],u);
				if(comps[ci][cii]==-1) {
					std::swap(comps[ci][cii],comps[ci].back());
					comps[ci].pop_back();
				} else {
					compof[comps[ci][cii]] = ci;
				}
				found=1;
				break;
			}
			if(!found && u!=-1) {
				comps[ci].push_back(u);
				compof[u] = ci;
				u=-1;
			}
			// if(u!=-1) compof[u]=-1;
		}
		assert(u==-1);
		for(int ci=0; ci<(int)comps.size(); ++ci) if(comps[ci].empty()) {
			comps[ci].swap(comps.back());
			comps.pop_back();
			for(auto u: comps[ci]) compof[u]=ci;
			--ci;
		}
	}
	return -imp;
}

// some stuff that can be optimized:
// - don't recompute the whole g every time
// - keep comps as an nxc matrix so you don't need to iterate colors

/// tries to move some nodes between components,
/// modifies comps and returns how much the score improved
int local_search(const problem& p, std::vector<std::vector<int>>& comps) {
	std::vector<int> compof(p.n(),-1);
	std::vector<std::vector<int>> og(p.n());
	for(int ci=0; ci<(int)comps.size(); ++ci)
		for(const auto u: comps[ci])
			compof[u]=ci;
	for(const auto [u,v]: p.edges) {
		og[u].push_back(v);
		og[v].push_back(u);
	}
	int C = *std::max_element(p.color.begin(), p.color.end())+1;
	int imp=0;
	bool good0=0;
	do {
		good0=0;
		for(int c=0; c<C; ++c) {
			bool good1=0;
			do{
				good1=0;
				int impc = local_search_color(p, og, compof, c, comps);
				if(impc) {
					good1=1;
					good0=1;
					imp+=impc;
				}
			} while(good1);
		}
	} while(good0);
	return imp;
}

