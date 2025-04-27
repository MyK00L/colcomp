#pragma once

#include <cmath>
#include <vector>
#include <Highs.h>
#include "mgraph.cpp"
#include "common.cpp"
#include "pricing.hpp"

/*
n_0 n_1 n_2 n_3...n_i
forall i: sum_j x[i,j]=n_i
forall j1,j2 max_i min(x[i,j1],x[i,j2])=1
max z = sum_j /\_i(x[i,j])

doable with SAT, I guess
but that's not the problem with what I'm doing anyway...
*/
std::vector<column> exact_pricing(const pricing_problem& p) {
    Highs h;
    // h.setOptionValue("seed",42); // TODO seed
    h.setOptionValue("output_flag", false);
    //h.setMatrixFormat(MatrixFormat::kColwise);
    h.changeObjectiveSense(ObjSense::kMaximize);
    int coli=0;
    for(size_t i=0; (int)i<p.g.n(); ++i) {
        h.addVar(0,1);
        h.changeColIntegrality(coli, HighsVarType::kInteger);
        // h.changeColCost(coli, -duals[i]);
        ++coli;
    }
    for(int i=0; i<int(p.g.n()); ++i) {
        for(int j=i+1; j<int(p.g.n()); ++j) {
            if(!compatible(p.g.c[i],p.g.c[j])) {
                int ind[2] = {i,j};
                double val[2] = {1,1};
                h.addRow(0,1,2,ind,val);
            }
        }
    }
    for(int i=0; i<int(p.g.n()); ++i) {
        for(auto [j,w]: p.g.g[i]) if(j>i && w>0) {
            h.addVar(0,1);
            //h.changeColIntegrality(coli, HighsVarType::kImplicitInteger);
            int nz[2] = {i,coli};
            double val[2] = {-1,1};
            h.addRow(-1,0,2,nz,val);
            nz[0]=j;
            h.addRow(-1,0,2,nz,val);
            h.changeColCost(coli,w);
            ++coli;
        }
    }
    // TODO integer cuts
    

    for(size_t i=0; i<p.cost.size(); ++i)
        h.changeColCost((int)i, -p.cost[i]);
    h.run();
    if(h.getObjectiveValue()<=0) return {};
    double obj = h.getObjectiveValue();
    std::vector<int> nz;
    std::vector<double> cv = h.getSolution().col_value;
    for(size_t i=0; i<p.cost.size(); ++i) {
        if(cv[i]>0.8) nz.push_back(i);
    }
    return {column{obj,nz}};

}
