#pragma once

#include <random>
#include <vector>
#include <array>
#include <utility>
#include "mgraph.cpp"
#include "heuristic_merge.cpp"

struct ShkDsu {
	ShkDsu(const std::vector<std::vector<int> >& colors):ufp(colors.size()) {
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
	bool try_merge(int u, int v) {
		u = find(u);
		v = find(v);
		if(cols[u].has_intersection(cols[v])) return 0;
		cols[u]|=cols[v];
		ufp[v]=u;
		return 1;
	}
};
std::pair<int,std::vector<int> > shk_score(const mgraph& mg, const std::vector<std::array<int,2> >& edg) {
	ShkDsu dsu(mg.c);
	std::vector<int> taken;
	for(int i=0; i<(int)edg.size(); ++i) if(dsu.try_merge(edg[i][0],edg[i][1])) taken.push_back(i);
	std::vector<std::vector<int> > sol(mg.n());
	for(int i=0; i<mg.n(); ++i) sol[dsu.find(i)].push_back(i);
	int score=0;
	for(auto &si: sol) if(!si.empty()) score+=mg.value_for_col(si);
	return {score, taken};
}
std::pair<int, std::vector<std::vector<int> > > single_heuristic_karger(mgraph mg, int nit) {
	std::vector<std::array<int,2> > edg;
	{
		std::vector<std::array<int,3> > tmp;
		for(int u=0; u<mg.n(); ++u) for(auto &[v,w]:mg.g[u]) if(v>u) {
			for(int _=0; _<w; ++_) tmp.push_back({cut_cost(mg,u,v)-merge_cost(mg,u,v),u,v});
		}
		sort(tmp.rbegin(),tmp.rend());
		edg.reserve(tmp.size());
		for(auto &[_,u,v]:tmp) edg.push_back({u,v});
	}
	Rng rng(4);
	auto [score,taken] = shk_score(mg,edg);

	for(int _=0; _<nit; ++_) {
		int i = taken[rng()%taken.size()];
		int j = rng()%edg.size();
		swap(edg[i],edg[j]);
		auto [nscore, ntaken] = shk_score(mg,edg);
		if(nscore<score) swap(edg[i],edg[j]);
		else {
			score=nscore;
			ntaken=taken;
		}
	}

	std::vector<std::vector<int> > sol(mg.n());
	{
		ShkDsu dsu(mg.c);
		std::vector<int> taken;
		for(int i=0; i<(int)edg.size(); ++i) if(dsu.try_merge(edg[i][0],edg[i][1])) taken.push_back(i);
		for(int i=0; i<mg.n(); ++i) for(auto oi: mg.groups[i]) sol[dsu.find(i)].push_back(oi);
		sol.erase(std::remove_if(sol.begin(),sol.end(),[](const std::vector<int>& v){return v.empty();}),sol.end());
	}
	
	return {score,sol};
}

std::vector<std::vector<int> > heuristic_karger(const problem& p) {
	mgraph mg(p);
	auto ccs = mg.split_components();
	int score=0;
	std::vector<std::vector<int> > ans;
	for(auto &cc: ccs) {
		auto [sc,sol] = single_heuristic_karger(cc,2000);
		score+=sc;
		std::move(sol.begin(),sol.end(),std::back_inserter(ans));
	}
	std::cerr<<"karger heuristic: "<<score<<std::endl;
	return ans;
}

