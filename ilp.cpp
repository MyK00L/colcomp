#include <Highs.h>
#include "common.cpp"

void ilp(problem p) {
    int n = (int)p.n();
    int m = (int)p.edges.size();
    int nvar_part = n+m;
    Highs highs;
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("seed",42); // TODO seed
    highs.changeObjectiveSense(ObjSense::kMaximize);
    // add variables
    for(int parti=0; parti<n; ++parti) {
        for(int ni=0; ni<n; ++ni) {
            highs.addVar(0,1);
            highs.changeColIntegrality(parti*nvar_part+ni, HighsVarType::kInteger);
        }
        for(int ei=0; ei<m; ++ei) {
            highs.addVar(0,1);
            highs.changeColCost(nvar_part*parti+n+ei,1);
        }
    }
    // add constraints
    for(int parti=0; parti<n; ++parti) {
        for(int i=0; i<n; ++i) for(int j=i+1; j<n; ++j) if(p.color[i]==p.color[j]) { // color constraints
            int ind[2] = {parti*nvar_part+i,parti*nvar_part+j};
            double val[2] = {1,1};
            highs.addRow(0,1,2,ind,val);
        }
        for(int ei=0; ei<m; ++ei) { // edge objective constraints
            int i = p.edges[ei][0];
            int j = p.edges[ei][1];
            int nz[2] = {parti*nvar_part+i,nvar_part*parti+n+ei};
            double val[2] = {-1,1};
            highs.addRow(-1,0,2,nz,val);
            nz[0]=parti*nvar_part+j;
            highs.addRow(-1,0,2,nz,val);
        }
    }
    for(int i=0; i<n; ++i) { // partition constraints
        vector<int> ind(n);
        vector<double> val(n,1);
        for(int parti=0; parti<n; ++parti) ind[parti]=parti*nvar_part+i;
        highs.addRow(0,1,n,ind.data(),val.data());
    }
    highs.run();
    std::cerr<<"ilp solution: "<<highs.getObjectiveValue()<<std::endl;
}