#pragma once

#include "common.cpp"
#include "mgraph.cpp"
#include <algorithm>

// FIXME
/*
struct beamsearchstate {
    beamsearchstate(mgraph _mg): mg(_mg), dual_bound(0) {
        std::vector<int> cfreq(1,0);
        dual_bound=0;
        for(int i=0; i<mg.n(); ++i) {
            for(auto &col: mg.c[i]) {
                if(col>=(int)cfreq.size()) cfreq.resize(col+1,0);
                cfreq[col]+=1;
            }
            for(auto &[j,w]:mg.g[i]) if(j>i) dual_bound+=w;
        }
        dual_bound-=*std::max_element(cfreq.begin(),cfreq.end());
    }
    mgraph mg;
    int dual_bound;
    std::array<beamsearchstate,2> split() const {
        mgraph mg1(mg),mg2(mg);
        int u=0; // TODO: better choice for split
        int v=mg.g[u][0][0];
        mg1.merge_edge(u,v);
        mg2.remove_edge(u,v);
        return {beamsearchstate(mg1), beamsearchstate(mg2)};
    }
    bool over() const {
        return mg.n()==0;
    }
};

std::vector<std::vector<int> > beamsearch(problem p, size_t beamsize=128) {
    int primal_bound = 0;
    mgraph best(p);
    std::vector<beamsearchstate> beam;
    beam.emplace_back(mgraph(p));
    while(!beam.empty()) {
        // std::cerr<<beam.size()<<std::endl;
        std::vector<beamsearchstate> newbeam;
        for(auto s = std::move_iterator(beam.begin()); s!=std::move_iterator(beam.end()); ++s) {
            if(s->dual_bound>primal_bound) {
                if(s->over()) {
                    best=s->mg;
                    primal_bound = s->dual_bound;
                } else {
                    for(auto s2: s->split()) newbeam.push_back(s2);
                }
            }
        }
        sort(newbeam.begin(), newbeam.end(), [](const auto &a, const auto &b){return a.dual_bound>b.dual_bound;});
        if(newbeam.size()>beamsize) newbeam.erase(newbeam.begin()+beamsize, newbeam.end());
        beam.swap(newbeam);
    }
    std::cerr<<"beamsearch: "<<primal_bound<<std::endl;
    return best.final_groups;
}
*/
