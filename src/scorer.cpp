#pragma once

#include "common.cpp"
#include <algorithm>
#include <iostream>

/// returns true if the solution is a valid solution for the problem
bool check(const problem &p, const std::vector<std::vector<int>> &sol) {
	std::vector<bool> vis(p.n(), 0);
	int C = *std::max_element(p.color.begin(), p.color.end())+1;
	for(auto &part : sol) {
		std::vector<bool> visc(C, 0);
		for(auto &u : part) {
			if(vis[u]) {
				std::cerr<<"invalid solution: same node present twice"<<std::endl;
				return 0;
			}
			if(visc[p.color[u]]) {
				std::cerr<<"invalid solution: same color present twice in the same partition"<<std::endl;
				return 0;
			}
			vis[u] = 1;
			visc[p.color[u]] = 1;
		}
	}
	/*
	for(int i=0; i<(int)p.n(); ++i) {
		if(vis[i]==0) {
			std::cerr<<"invalid solution: not all nodes are present"<<std::endl;
			return 0;
		}
	}
	*/
	return 1;
}

int score(const problem &p, const std::vector<std::vector<int>> &sol) {
	if(!check(p, sol)) return -1;
	std::vector<int> part(p.n(), -1);
	for(int pi = 0; pi<(int)sol.size(); ++pi) {
		for(auto &u : sol[pi]) part[u] = pi;
	}
	int ans = 0;
	for(auto &[u, v] : p.edges) if(part[u]!=-1) { ans += part[u]==part[v]; }
	return ans;
}
