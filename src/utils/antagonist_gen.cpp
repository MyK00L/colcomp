#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>

using namespace std;

/*
based on the following graph:
v 0 0
v 1 0
v 2 0
v 3 1
v 4 1
v 5 2
v 6 2
v 7 3
e 0 3
e 1 3
e 2 4
e 2 6
e 2 7
e 3 5
e 3 7
e 4 6
e 5 7
*/
const int cbad[8] = {0,0,0,1,1,2,2,3};
const int ebad[9][2] = {
	{0,3},
	{1,3},
	{2,4},
	{2,6},
	{2,7},
	{3,5},
	{3,7},
	{4,6},
	{5,7}
};
/// adds the above graph, with random colors
void add_badgraph(mt19937& rng, vector<int>& color, vector<array<int,2>>& edg, const int n, const int c) {
	vector<bool> visc(c);
	array<int,4> cs;
	for(int i=0; i<4; ++i) {
		cs[i]=rng()%(c-i);
		while(visc[cs[i]]) ++cs[i];
		visc[cs[i]]=1;
	}
	int u0=color.size();
	for(int i=0; i<8; ++i) color.push_back(cs[cbad[i]]);
	for(int ei=0; ei<9; ++ei) {
		if(u0+ebad[ei][1]<n)
			edg.push_back({ebad[ei][0]+u0,ebad[ei][1]+u0});
	}
}

int main(int argv, char **args) {
	if(argv != 5) {
		cerr << "usage: " << string(args[0]) << " seed n c d\n";
		return 1;
	}
	int seed = atoi(args[1]);
	int n = atoi(args[2]);
	int c = atoi(args[3]);
	double density = atof(args[4]);
	mt19937 rng(seed);
	vector<int> color;
	vector<array<int, 2>> edg;

	// add the small bad graphs
	while(int(color.size())<n) add_badgraph(rng, color, edg, n, c);
	// connect them in a tree randomly so as to have a single cc
	for(int i=8; i<n; i+=8) {
		int ui = i+(rng()%min(n-i,8));
		int vi = rng()%i;
		while(color[ui]==color[vi]) vi=(vi+1)%i;
		edg.push_back({vi,ui});
	}

	set<array<int,2>> present_edg;
	for(auto &e: edg) {
		present_edg.insert(e);
		array<int,2> e1 = {e[1],e[0]};
		present_edg.insert(e1);
	}
	vector<array<int, 2>> extra_edgs;
	for(int i=0; i<n; ++i)
		for(int j=i+1; j<n; ++j)
			if(color[i]!=color[j] && !present_edg.count({i,j}))
				extra_edgs.push_back({i,j});
	shuffle(extra_edgs.begin(),extra_edgs.end(),rng);
	int em = min(int(ceil(extra_edgs.size() * density)), int(extra_edgs.size()));
	for(int i=0; i<em; ++i) edg.push_back(extra_edgs[i]);

	int m = edg.size();
	std::cerr<<"m="<<m<<std::endl;
	for(int i = 0; i < n; ++i) { cout << "v " << i << ' ' << color[i] << '\n'; }
	for(int ei = 0; ei < m; ++ei) { cout << "e " << edg[ei][0] << ' ' << edg[ei][1] << '\n'; }
	return 0;
}
