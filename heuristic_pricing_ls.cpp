#pragma once
#include <vector>
#include <array>

double swtch(std::vector<bool>& ans, const std::vector<double>& duals, const std::vector<std::vector<std::array<int,2> > >& g, const int i) {
    double res;
    if(ans[i]) {
        res=duals[i];
        for(const auto &[j,w]: g[i]) if(ans[j]) res-=w;
    } else {
        res=-duals[i];
        for(const auto &[j,w]: g[i]) if(ans[j]) res+=w;
    }
    ans[i].flip();
    return res;
}

double local_search_pricing(std::vector<bool>& ans, const std::vector<double>& duals, const std::vector<std::vector<std::array<int,2> > >& g) {
    double res=0;
    bool up;
    do {
        up=0;
        for(int i=0; i<(int)ans.size(); ++i) {
            double p = swtch(ans, duals, g, i);
            if(p>0) {
                up=1;
                res+=p;
            } else {
                swtch(ans,duals,g,i);
            }
        }
    } while(up);
    return res;
}
