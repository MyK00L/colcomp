#include <array>
#include <ctime>
#include <iostream>

#include "common.cpp"
#include "heuristic_beamsearch.cpp"
#include "bnb.cpp"
#include "bnp.cpp"
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

int main() {
	cin.tie(0);
	ios_base::sync_with_stdio(0);
	
	problem p = input();
	// Master m(p);
	// m.run();
	clock_t t0 = clock();
	//ilp(p);
	clock_t t1 = clock();
	bnp_main(p);
	clock_t t2 = clock();
	std::cerr<<"ilp time: "<<double(t1-t0)/CLOCKS_PER_SEC<<'\n';
	std::cerr<<"bnp time: "<<double(t2-t1)/CLOCKS_PER_SEC<<'\n';
	
	return 0;
}

// heuristic idea: keep cutting in the min-cut between two nodes with the same color
