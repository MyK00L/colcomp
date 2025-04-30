#pragma once

#include <random>
#include <vector>
#include <array>
#include "mgraph.cpp"

struct ShkDsu {
	ShkDsu(std::vector<std::vector<int> > colors):ufp(colors.size()) {
		int C=0;
		for(auto&cs:colors) for(auto &c:cs) if(c+1>C) C=c+1;
		cols = std::vector<bitvec>(colors.size(), bitvec(C));
		std::iota(ufp.begin(),ufp.end(),0);
		for(int i=0; i<(int)cols.size(); ++i) {
			for(auto &c:colors[i]) cols[i].set1(c);
		}
	}
	std::vector<int> ufp; // parent in ufs
	std::vector<bitvec> cols; // colors
	int find(int u) {
		return u==ufp[u]?u:ufp[u]=find(ufp[u]);
	}
	void try_merge(int u, int v) {
		u = find(u);
		v = find(v);
		if(cols[u].has_intersection(cols[v])) return;
		cols[u]|=cols[v];
		ufp[v]=u;
	}
};
int single_heuristic_karger(mgraph mg, int nit) {
	std::vector<std::array<int,2> > edg;
	for(int u=0; u<mg.n(); ++u) for(auto &[v,w]:mg.g[u]) if(v>u) {
		for(int _=0; _<w; ++_) edg.push_back({u,v});
	}
	std::mt19937 rng(4);
	const ShkDsu dsu0(mg.c);
	int best_score=0;
	while(nit--) {
		ShkDsu dsu = dsu0;
		std::shuffle(edg.begin(), edg.end(), rng);
		for(auto &[u,v]: edg) {
			dsu.try_merge(u,v);
		}
		std::vector<std::vector<int> > sol(mg.n());
		for(int i=0; i<mg.n(); ++i) sol[dsu.find(i)].push_back(i);
		int score=0;
		for(auto &si: sol) if(!si.empty()) {
			score+=mg.value_for_col(si);
		}
		if(score>best_score) best_score=score;
	}
	return best_score;
}

void heuristic_karger(const problem& p) {
	mgraph mg(p);
	auto ccs = mg.split_components();
	int score=0;
	for(auto &cc: ccs) {
		score+=single_heuristic_karger(cc,10000);
	}
	std::cerr<<"karger heuristic: "<<score<<std::endl;
}

