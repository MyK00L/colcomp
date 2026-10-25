#pragma once

#include "mgraph.cpp"
#include <array>
#include <iostream>
#include <queue>
#include <vector>

struct pricing_problem {
	mgraph g;
	std::vector<double> cost;
	bool is_node_useless(int u, const std::vector<bool> &nu, const std::vector<std::vector<bool>> &eu) {
		std::vector<int> ma(g.C(), 0);
		for(auto &[v, w] : g.g[u])
			if(!nu[v] && !eu[u][v]) ma[g.c[v][0]] = std::max(ma[g.c[v][0]], w);
		double x = -cost[u];
		for(auto &y : ma) x += y;
		return x<0;
	}
	bool is_edge_useless(int u, int v, const std::vector<bool> &nu, const std::vector<std::vector<bool>> &eu) {
		std::vector<int> ma(g.C(), 0);
		{
			auto ui = g.g[u].begin();
			auto vi = g.g[v].begin();
			while(ui!=g.g[u].end() || vi!=g.g[v].end()) {
				if(vi!=g.g[v].end() && (ui==g.g[u].end() || (*ui)[0]>(*vi)[0])) {
					int j = (*vi)[0];
					if(!nu[j] && !eu[v][j]) {
						int w = (*vi)[1];
						ma[g.c[j][0]] = std::max(ma[g.c[j][0]], w);
					}
					++vi;
				} else if(ui!=g.g[u].end() && (vi==g.g[v].end() || (*vi)[0]>(*ui)[0])) {
					int j = (*ui)[0];
					if(!nu[j] && !eu[u][j]) {
						int w = (*ui)[1];
						ma[g.c[j][0]] = std::max(ma[g.c[j][0]], w);
					}
					++ui;
				} else {
					int j = (*ui)[0];
					if(!nu[j] && !eu[u][j] && !eu[v][j]) {
						int w = (*ui)[1]+(*vi)[1];
						ma[g.c[j][0]] = std::max(ma[g.c[j][0]], w);
					}
					++ui;
					++vi;
				}
			}
		}
		int uvw = 0;
		for(auto &[j, w] : g.g[u])
			if(j==v) uvw = w;
		for(auto &cu : g.c[u]) ma[cu] = 0;
		for(auto &cv : g.c[v]) ma[cv] = 0;
		double x = -cost[u]-cost[v]+uvw;
		for(auto &y : ma) x += y;
		return x<0;
	}
	// WRONG, to fix or delete (eg on BBS11037C5.graph)
	void preprocess() {
		std::vector<bool> nu(g.n(), 0);										   // 1 if node is useless
		std::vector<std::vector<bool>> eu(g.n(), std::vector<bool>(g.n(), 0)); // 1 if edge is useless
		std::vector<bool> inq(g.n(), 1);
		std::queue<int> q;
		for(int i = 0; i<g.n(); ++i) q.push(i);
		bool upd = 1;
		while(upd) {
			while(!q.empty()) {
				int u = q.front();
				q.pop();
				if(is_node_useless(u, nu, eu)) {
					nu[u] = 1;
					std::vector<std::array<int, 2>> tredg;
					for(auto &[v, _w] : g.g[u]) {
						eu[u][v] = eu[v][u] = 1;
						tredg.push_back({u, v});
						if(!nu[v] && !inq[v]) {
							inq[v] = 1;
							q.push(v);
						}
					}
				}
				inq[u] = 0;
			}
			upd = 0;
			for(int u = 0; u<g.n() && !upd; ++u)
				if(!nu[u])
					for(auto &[v, w] : g.g[u])
						if(!nu[v] && !eu[u][v]) {
							if(is_edge_useless(u, v, nu, eu)) {
								upd = 1;
								eu[u][v] = eu[v][u] = 1;
								inq[u] = 1;
								inq[v] = 1;
								q.push(u);
								q.push(v);
								break;
							}
						}
		}
		for(int i = 0; i<g.n(); ++i) {
			for(auto &[j, w] : g.g[i])
				if(eu[i][j]) w = -g.orig_n*g.orig_n-1;
			if(nu[i]) cost[i] = g.orig_n+1;
		}
	}
	pricing_problem(mgraph _g, std::vector<double> duals) : g(_g), cost(duals) { preprocess(); }
};

inline void checkcolumn(const mcolumn &col, const pricing_problem &pp) {
	std::vector<bool> visc(pp.g.C(), 0);
	for(const auto &u : col.nodes) {
		for(const auto &c : pp.g.c[u]) {
			assert(!visc[c]);
			visc[c] = 1;
		}
	}
}
inline void debugcolumn(const mcolumn &col, const pricing_problem &pp) {
	double val = col.value;
	for(const auto &u : col.nodes) val -= pp.cost[u];
	std::cerr<<"colval: "<<val<<std::endl;
	std::vector<bool> inc(pp.g.n());
	for(const auto &u : col.nodes) {
		inc[u] = 1;
		std::cerr<<u<<"[-"<<pp.cost[u]<<"](";
		for(const auto &c : pp.g.c[u]) std::cerr<<c<<",";
		std::cerr<<") ";
	}
	std::cerr<<std::endl;
	for(const auto &u : col.nodes) {
		for(const auto &v : col.nodes)
			if(v>u) { std::cerr<<u<<"-"<<v<<": "<<pp.g.edge_val(u, v)<<std::endl; }
	}
}
