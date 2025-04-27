#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include "heuristic_pricing_ls.cpp"
#include "pricing.hpp"

std::vector<mcolumn> heuristic_pricing_annealing(const pricing_problem& pp) {
    int C=pp.g.C();
    const std::vector<double>& duals = pp.cost;
    const int n = pp.g.n();
    std::vector<bool> best_ans(n,0);
    double best_val = 0;
    std::vector<bool> ans = best_ans;
    double val = 0;
    std::vector<int> by_color(C,-1);
    double temp = 64;
    Rng rng(4);
    while(temp>=0.1) {
        int i = int(rng()%n);
        double new_val = val;
        std::vector<int> to_switch;
        to_switch.push_back(i);
        if(!ans[i]) {
            for(auto &ci: pp.g.c[i]) if(by_color[ci]!=-1) to_switch.push_back(by_color[ci]);
        }
        for(auto &j: to_switch) new_val+=swtch(ans,duals,pp.g.g,j);
        if(new_val>val || rng.f01()>std::exp((new_val-val)/temp)) {
            val = new_val;
            if(ans[i]) for(auto &ci: pp.g.c[i]) by_color[ci]=i;
            else for(auto &ci: pp.g.c[i]) by_color[ci]=-1;
            if(val>best_val) {
                best_val=val;
                best_ans=ans;
            }
        } else {
            for(auto &j: to_switch) swtch(ans,duals,pp.g.g,j);
        }
        temp*=0.999;
    }
    best_val+=local_search_pricing(best_ans,duals,pp.g.g);
    if(best_val<=0.0001) return {};
    std::vector<int> nz;
    for(int i=0; i<n; ++i) if(best_ans[i]) nz.push_back(i);
    return {mcolumn{pp.g.value_for_col(nz),nz}};
}
