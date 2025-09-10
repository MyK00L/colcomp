#include <array>
#include <ctime>
#include <iostream>

#include "common.cpp"
#include "heuristic/heuristic_beamsearch.cpp"
#include "heuristic/heuristic_karger.cpp"
#include "heuristic/heuristic_merge.cpp"
#include "bnb.cpp"
#include "bnp/bnp.cpp"
#include "ilp.cpp"

using namespace std;

/*
 *	MCC: minimize number of components min cn
 *	MEC: maximize transitive closure sum(cs*(cs-1))
 *	MOP: minimize number of removed edges sum(cut(c))
 */

/*
 * 2 colors -> bipartite matching
 *
 */

problem input() {
	vector<uint32_t> col;
	vector<array<uint32_t,2> > edg;
	while(!cin.eof()) {
		char c=0;
		cin>>c;
		if(c=='v') {
			uint32_t a,b;
			cin>>a>>b;
			if(a+1>col.size()) col.resize(a+1);
			col[a]=b;
		} else if(c=='e') {
			uint32_t a,b;
			cin>>a>>b;
			edg.push_back({a,b});
		} else {
			break;
		}
	}
	return problem{col,edg};
}

int ccs(mgraph mg) {
	vector<bool> vis(mg.n(),0);
	int res=0;
	for(int i=0; i<mg.n(); ++i) if(!vis[i]) {
		res+=1;
		queue<int> q;
		q.push(i);
		vis[i]=1;
		while(!q.empty()) {
			int u = q.front();
			q.pop();
			for(auto &[v,_w]:mg.g[u]) if(!vis[v]) {
				vis[v]=1;
				q.push(v);
			}
		}
	}
	return res;
}

int main() {
	cin.tie(0);
	ios_base::sync_with_stdio(0);

	problem p = input();
	std::cerr<<"|E|="<<p.edges.size()<<std::endl;
	std::cerr<<"ccs="<<ccs(mgraph(p))<<std::endl;
	// Master m(p);
	// m.run();
	clock_t t0 = clock();
	// ilp(p);
	heuristic_karger(p);
	heuristic_merge(p);
	clock_t t1 = clock();
	bnp_main(p);
	clock_t t2 = clock();
	std::cerr<<"ilp time: "<<double(t1-t0)/CLOCKS_PER_SEC<<'\n';
	std::cerr<<"bnp time: "<<double(t2-t1)/CLOCKS_PER_SEC<<'\n';

	std::cerr<<"nnodes: "<<nnodes<<std::endl;
	
	return 0;
}

// heuristic idea: keep cutting in the min-cut between two nodes with the same color
