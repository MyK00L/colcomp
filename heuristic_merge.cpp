#include <vector>
#include <iostream>

#include "mgraph.cpp"
#include "common.cpp"

int merge_cost(const mgraph& mg, int u, int v) {
	int ans = 0;
	for(auto &[j,w]: mg.g[u]) if(!compatible(mg.c[v],mg.c[j])) ans += w;
	for(auto &[j,w]: mg.g[v]) if(!compatible(mg.c[u],mg.c[j])) ans += w;
	return ans;
}
int cut_cost(const mgraph& mg, int u, int v) {
	int ans = 3*mg.edge_val(u,v);
	auto ui = mg.g[u].begin();
	auto ue = mg.g[u].end();
	auto vi = mg.g[v].begin();
	auto ve = mg.g[v].end();
	while(ui!=ue && vi!=ve) {
		if(ui->at(0) < vi->at(0)) ++ui;
		else if(vi->at(0) < ui->at(0)) ++vi;
		else {
			// ans += std::min(ui->at(1),vi->at(1));
			ans += ui->at(1) + vi->at(1);
			++ui;
			++vi;
		}
	}
	return ans;
}

int heuristic_merge_comp(mgraph mg) {
	int score=0;
	while(1) {
		int bu=-1,bv=-1,best=std::numeric_limits<int>::min();
		for(int u=0; u<mg.n(); ++u) {
			for(auto [v,_w]: mg.g[u]) if(v>u) {
				int t = cut_cost(mg,u,v) - merge_cost(mg,u,v);
				if(t>best) {
					best=t;
					bu=u;
					bv=v;
				}
			}
		}
		if(bu==-1) break;
		score += mg.merge_edge(bu,bv);
	}
	return score;
}

void heuristic_merge(const problem& p) {
	mgraph mg(p);
	auto ccs = mg.split_components();
	int score=0;
	for(auto cc:ccs) {
		score+=heuristic_merge_comp(cc);
	}
	std::cerr<<"heuristic merge: "<<score<<std::endl;
}

