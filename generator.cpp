#include <cstdlib>
#include <random>
#include <iostream>
#include <algorithm>
#include <array>

using namespace std;

int main(int argv, char** args) {
    if(argv!=5) {
        cerr<<"usage: "<<string(args[0])<<" seed n c density\n";
        return 1;
    }
    int seed = atoi(args[1]);
    int n = atoi(args[2]);
    int c = atoi(args[3]);
    double density = atof(args[4]);
    mt19937 rng(seed);
    vector<int> color(n);
    for(auto &i:color) i = rng()%c;
    sort(color.begin(), color.end());
    vector<array<int,2> > edg;
    for(int i=0; i<n; ++i) for(int j=i+1; j<n; ++j) if(color[i]!=color[j]) {
        edg.push_back({i,j});
    }
    shuffle(edg.begin(),edg.end(),rng);
    int m = min(int(ceil(edg.size()*density)),int(edg.size()));
    for(int i=0; i<n; ++i) {
        cout<<"v "<<i<<' '<<color[i]<<'\n';
    }
    for(int ei=0; ei<m; ++ei) {
        cout<<"e "<<edg[ei][0]<<' '<<edg[ei][1]<<'\n';
    }
    
    return 0;
}
