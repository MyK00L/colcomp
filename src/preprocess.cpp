#pragma once

#include "common.cpp"
#include "heuristic/heuristic_karger.cpp"
#include "heuristic/heuristic_merge.cpp"
#include <algorithm>
#include <iostream>
#include <queue>
#include <set>

// could be done faster with fast all-pairs minimum cut
/// given a subgraph, compute some node-node cuts
/// og is the original graph
/// comp the list of nodes of the subgraph
/// octs is the list of node pairs to compute cuts for (deduped inside)
/// O(|octs|*nm+nc) (where n=|comp| and m is the number of edges in it)
std::vector<int> compute_component_cuts(const std::vector<std::vector<int>> &og, std::vector<int> comp,
										const std::vector<std::array<int, 2>> &octs) {
	struct Edge {
		int v;
		int w;
		int rwj;
	};
	std::vector<std::vector<Edge>> g(comp.size()); // compressed graph
	sort(comp.begin(), comp.end());
	for(int ui = 0; ui<(int)comp.size(); ++ui) {
		const int u = comp[ui];
		for(auto &v : og[u])
			if(v>u) {
				auto it = std::lower_bound(comp.begin(), comp.end(), v);
				if(it!=comp.end() && *it==v) {
					const int vi = it-comp.begin();
					g[ui].push_back(Edge{vi, 1, int(g[vi].size())});
					g[vi].push_back(Edge{ui, 1, int(g[ui].size()-1)});
				}
			}
	}
	std::vector<std::array<int, 3>> cts;
	for(const auto &[uo, vo] : octs) {
		auto it = std::lower_bound(comp.begin(), comp.end(), uo);
		assert(*it==uo);
		int u = it-comp.begin();
		it = std::lower_bound(comp.begin(), comp.end(), vo);
		if(*it!=vo) {
			std::cerr<<"comp: ";
			for(auto x : comp) std::cerr<<x<<' ';
			std::cerr<<"\nocts: ";
			for(auto x : octs) std::cerr<<x[0]<<','<<x[1]<<' ';
			std::cerr<<std::endl;
		}
		assert(*it==vo);
		int v = it-comp.begin();
		if(u>v) std::swap(u, v);
		cts.push_back({u, v, -1});
	}
	sort(cts.begin(), cts.end());
	cts.erase(unique(cts.begin(), cts.end()), cts.end());
	std::vector<int> ff(g.size(), -1);
	std::vector<std::array<int *, 2>> fe(g.size());
	for(auto &[u0, u1, flow] : cts) {
		flow = 0;
		bool found = 1;
		while(found) {
			found = 0;
			std::queue<int> q;
			std::fill(ff.begin(), ff.end(), -1);
			q.push(u0);
			ff[u0] = u0;
			while(!q.empty()) {
				int u = q.front();
				q.pop();
				if(u==u1) {
					while(u!=ff[u]) {
						assert(*fe[u][0]>0);
						*fe[u][0] -= 1;
						*fe[u][1] += 1;
						u = ff[u];
					}
					found = 1;
					flow += 1;
					break;
				}
				for(auto &e : g[u])
					if(e.w>0 && ff[e.v]==-1) {
						ff[e.v] = u;
						fe[e.v] = {&e.w, &g[e.v][e.rwj].w};
						q.push(e.v);
					}
			}
		}
		// reset flow
		for(auto &el : g)
			for(auto &e : el) e.w = 1;
	}
	for(auto &[u0, u1, flow] : cts) {
		u0 = comp[u0];
		u1 = comp[u1];
		if(u0>u1) std::swap(u0, u1);
	}
	std::vector<int> ans;
	ans.reserve(octs.size());
	for(auto [ou, ov] : octs) {
		if(ou>ov) std::swap(ou, ov);
		std::array<int, 3> hi = {ou, ov, 0};
		auto it = std::lower_bound(cts.begin(), cts.end(), hi);
		assert(it->at(0)==ou && it->at(1)==ov);
		ans.push_back(it->at(2));
	}
	return ans;
}

// should be possible (probably hard and not worth it) to reduce the c^3 with fast all-pairs minimum cut
// should be easy to reduce the m^2 if needed
/// Try removing as many edges as you can from those not in comps using the rule
/// Complexity is O(sum_{comp}{|comp|^4}+m^2) = O(nc^3+m^2) (but it should be faster that what it sounds like)
int try_with_sol(problem &p, std::vector<std::vector<int>> &g, const std::vector<std::vector<int>> &comps,
				 const std::vector<bool> &ignore) {
	std::queue<int> q; // q of components to check
	std::vector<bool> inq(comps.size());
	std::set<std::array<int, 2>> removed; // set of removed edges
	std::vector<int> node_comp(p.n());	  // for each node, component it is in
	std::vector<std::vector<std::array<int, 3>>> comp_edges(
		comps.size()); // for each component, its adjacent edges with the weight of the cut needed
	const int C = *std::max_element(p.color.begin(), p.color.end())+1;

	auto compute_cuts_for = [&](int compi) {
		std::vector<int> col_node(C, -1);
		for(auto u : comps[compi]) col_node[p.color[u]] = u;
		std::vector<std::array<int, 2>> tocheck;
		for(auto [_, u, v] : comp_edges[compi]) {
			if(node_comp[u]!=compi) std::swap(u, v); // don't swap them in comp_edges
			if(col_node[p.color[v]]!=-1) tocheck.push_back({u, col_node[p.color[v]]});
		}
		std::vector<int> ws = compute_component_cuts(g, comps[compi], tocheck);
		int wi = 0;
		for(auto &[w, u, v] : comp_edges[compi]) {
			const int cv = node_comp[u]==compi ? col_node[p.color[v]] : col_node[p.color[u]];
			if(cv!=-1) w = ws[wi++];
			else w = -1;
		}
		sort(comp_edges[compi].begin(), comp_edges[compi].end());
	};

	for(int compi = 0; compi<(int)comps.size(); ++compi)
		if(!ignore[comps[compi][0]]) {
			std::vector<int> col_node(C, -1);
			for(auto u : comps[compi]) {
				col_node[p.color[u]] = u;
				node_comp[u] = compi;
			}
			for(auto u : comps[compi]) {
				for(auto v : g[u])
					if(node_comp[v]!=node_comp[u]) {
						if(u<v) comp_edges[compi].push_back({-1, u, v});
						else comp_edges[compi].push_back({-1, v, u});
					}
			}
			// only compute if |comp|>|edges|, otherwise you can't remove any edge anyway
			if(comps[compi].size()>comp_edges[compi].size()) {
				compute_cuts_for(compi);
				q.push(compi);
				inq[compi] = 1;
			}
		}

	// do the thing
	while(!q.empty()) {
		int compi = q.front();
		q.pop();
		bool good_before = comps[compi].size()>comp_edges[compi].size();
		comp_edges[compi].erase(
			std::remove_if(comp_edges[compi].begin(), comp_edges[compi].end(),
						   [&removed](std::array<int, 3> e) { return removed.count({e[1], e[2]}); }),
			comp_edges[compi].end());
		if(comps[compi].size()>comp_edges[compi].size()) {
			if(!good_before) compute_cuts_for(compi);
			while(comp_edges[compi].size() && comp_edges[compi].back()[0]>=(int)comp_edges[compi].size()) {
				std::array<int, 2> edg = {comp_edges[compi].back()[1], comp_edges[compi].back()[2]};
				removed.insert(edg);
				for(int u : edg)
					if(!inq[node_comp[u]] && !ignore[comps[node_comp[u]][0]]) {
						q.push(node_comp[u]);
						inq[node_comp[u]] = 1;
					}
				comp_edges[compi].pop_back();
			}
			inq[compi] = 0;
		}
	}

	// remove from p
	for(int ei = 0; ei<(int)p.edges.size(); ++ei) {
		auto [u, v] = p.edges[ei];
		if(u>v) std::swap(u, v);
		std::array<int, 2> edg = {int(u), int(v)};
		if(removed.count(edg)) {
			std::swap(p.edges[ei], p.edges.back());
			p.edges.pop_back();
			--ei;
		}
	}
	// remove from g
	for(int u = 0; u<int(g.size()); ++u) {
		for(int ei = 0; ei<int(g[u].size()); ++ei) {
			std::array<int, 2> edg = {u, g[u][ei]};
			if(edg[0]>edg[1]) std::swap(edg[0], edg[1]);
			if(removed.count(edg)) {
				std::swap(g[u][ei], g[u].back());
				g[u].pop_back();
				--ei;
			}
		}
	}
	return removed.size();
}

void update_finished(const problem &p, const std::vector<std::vector<int>> &g, std::vector<bool> &finished) {
	const int C = *std::max_element(p.color.begin(), p.color.end())+1;
	std::vector<bool> vis = finished;
	for(int u0 = 0; u0<int(g.size()); ++u0)
		if(!vis[u0]) {
			std::vector<bool> visc(C);
			bool good = 1;
			std::vector<int> q;
			int qi = 0;
			q.push_back(u0);
			vis[u0] = 1;
			while(qi!=(int)q.size()) {
				int u = q[qi++];
				if(visc[p.color[u]]) good = 0;
				visc[p.color[u]] = 1;
				for(auto &v : g[u]) {
					if(!vis[v]) {
						vis[v] = 1;
						q.push_back(v);
					}
				}
			}
			if(good)
				for(auto &u : q) finished[u] = 1;
		}
}

std::vector<std::array<int, 2>> can_remove(const problem &p, const std::vector<std::vector<int>> &g,
										   const std::vector<int> comp, const std::vector<int> &col_node,
										   const std::vector<bool> &incomp) {
	std::vector<std::array<int, 2>> tocheck;
	std::vector<std::array<int, 3>> edg;
	int nedg = 0;
	for(auto &u : comp) {
		for(auto &v : g[u])
			if(!incomp[v]) {
				++nedg;
				if(col_node[p.color[v]]!=-1) {
					edg.push_back({-1, u, v});
					tocheck.push_back({u, col_node[p.color[v]]});
				}
			}
	}
	const std::vector<int> ws = compute_component_cuts(g, comp, tocheck);
	for(int ei = 0; ei<(int)edg.size(); ++ei) edg[ei][0] = ws[ei];
	std::sort(edg.rbegin(), edg.rend());
	std::vector<std::array<int, 2>> ans;
	for(int i = 0; i<(int)edg.size() && edg[i][0]>=nedg-i; ++i) {
		int u = edg[i][1], v = edg[i][2];
		if(u>v) std::swap(u, v);
		ans.push_back({u, v});
	}
	return ans;
}
int try_with_greedy(problem &p, std::vector<std::vector<int>> &g, const std::vector<bool> &ignore) {
	const int C = *std::max_element(p.color.begin(), p.color.end())+1;
	std::vector<bool> incomp(p.n());
	std::set<std::array<int, 2>> removed;
	for(int u0 = 0; u0<int(p.n()); ++u0)
		if(!ignore[u0]) {
			int nedges = 0;
			std::vector<int> comp;
			std::set<int> neigh;
			std::vector<int> col_node(C, -1);
			comp.push_back(u0);
			incomp[u0] = 1;
			col_node[p.color[u0]] = u0;
			for(int v : g[u0]) {
				neigh.insert(v);
				nedges += 1;
			}
			for(;;) {
				int minn = p.n()*2;
				int minu = -1;
				for(const auto u : neigh)
					if(col_node[p.color[u]]==-1) {
						int unn = 0;
						for(const auto v : g[u]) {
							if(incomp[v]) unn -= 1;
							else unn += 1;
						}
						if(unn<minn) {
							minn = unn;
							minu = u;
						}
					}
				if(minu==-1) break;
				nedges += minn;
				col_node[p.color[minu]] = minu;
				comp.push_back(minu);
				incomp[minu] = 1;
				neigh.erase(minu);
				for(const auto v : g[minu])
					if(!incomp[v]) neigh.insert(v);

				if(nedges<(int)comp.size()) { // cut would be always<|edges| otherwise
					const auto cr = can_remove(p, g, comp, col_node, incomp);
					// remove from g, might need to optimize this
					for(const auto &[u, v] : cr) {
						if(u>v) removed.insert({v, u});
						else removed.insert({u, v});
						for(int ei = 0; ei<(int)g[u].size(); ++ei)
							if(g[u][ei]==v) {
								std::swap(g[u][ei], g[u].back());
								g[u].pop_back();
								break;
							}
						for(int ei = 0; ei<(int)g[v].size(); ++ei)
							if(g[v][ei]==u) {
								std::swap(g[v][ei], g[v].back());
								g[v].pop_back();
								break;
							}
					}
				}
			}
			for(const auto u : comp) incomp[u] = 0;
		}

	// remove from p
	for(int ei = 0; ei<(int)p.edges.size(); ++ei) {
		auto [u, v] = p.edges[ei];
		if(u>v) std::swap(u, v);
		std::array<int, 2> edg = {int(u), int(v)};
		if(removed.count(edg)) {
			std::swap(p.edges[ei], p.edges.back());
			p.edges.pop_back();
			--ei;
		}
	}
	return removed.size();
}

// TODO: mark and ignore those in finished components,
//       try with both heuristics, if fail try greedily constructed sets (like before change)
void preprocess(problem &p) {
	std::vector<bool> finished(p.n());
	std::vector<std::vector<int>> g(p.n());
	for(auto &[u, v] : p.edges) {
		g[u].push_back(v);
		g[v].push_back(u);
	}

	update_finished(p, g, finished);
	try_with_greedy(p, g, finished);
	update_finished(p, g, finished);
	try_with_sol(p, g, heuristic_merge(p), finished);
	int removed = 0;
	do {
		removed = 0;
		update_finished(p, g, finished);
		removed += try_with_sol(p, g, heuristic_karger(p, 1000), finished);
	} while(removed);
	update_finished(p, g, finished);
	removed = try_with_greedy(p, g, finished);
}
