#include <iostream>
#include <cstdint>
#include <vector>
#include <array>

using namespace std;

struct problem {
	std::vector<uint32_t> color;
	std::vector<std::array<uint32_t,2> > edges;
	size_t n() {
		return color.size();
	}
};


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

/*
param n := 6;
param m := 7;

set E := (1,2) (2,3) (3,4) (4,5) (5,6) (1,6) (1,3);

param col :=
1 1
2 1
3 2
4 3
5 2
6 3
;
*/
int main() {
    problem p = input();
    cout<<"param n := "<<p.n()<<";\n";
    cout<<"param m := "<<p.edges.size()<<";\n";
    cout<<"set E :=";
    for(auto [u,v]: p.edges) cout<<" ("<<u+1<<","<<v+1<<")";
    cout<<";\n";
    cout<<"param col :=\n";
    for(int i=0; i<(int)p.n(); ++i) cout<<i+1<<' '<<p.color[i]+1<<'\n';
    cout<<";\n";
    return 0;
}
