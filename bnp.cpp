#pragma once

#include <map>
#include <iostream>
#include "common.cpp"
#include "mgraph.cpp"
#include "master.hpp"
#include "master.cpp"

std::vector<mcolumn> columns_for(const std::vector<mcolumn>& cols, const mgraph& mg) {
    std::vector<mcolumn> ans;
    // single-element columns for feasibility
    for(int i=0; i<mg.n(); ++i) ans.push_back({0,{i}});

    std::vector<int> group_of_node(mg.orig_n,-1); // what group each node is in in the subproblem, -1 if finalized
    for(int i=0; i<(int)mg.groups.size(); ++i) for(auto &u: mg.groups[i]) group_of_node[u]=i;
    std::vector<int> ng(mg.groups.size(),0);
    for(const auto& [_v,col]: cols) {
        bool ok=1; // if column is compatible with subproblem
        for(auto &u:col) {
            if(group_of_node[u]<0) {
                ok=0;
                break;
            }
            ng[group_of_node[u]]+=1;
        }
        if(ok) for(auto &u: col) {
            if(ng[group_of_node[u]]!=(int)mg.groups[group_of_node[u]].size()) {
                ok=0;
                break;
            }
        }
        if(ok) { // add column
            int val=0;
            std::vector<int> mgcol;
            for(auto &uu: col) if(mg.groups[group_of_node[uu]][0]==uu) {
                int u = group_of_node[uu];
                mgcol.push_back(u);
                for(auto &[v,w]: mg.g[u]) if(v>u && ng[v]) {
                    val+=w;
                }
            }
            ans.push_back({val,mgcol});
        }
        // reset ng
        for(auto &u:col) {
            if(group_of_node[u]<0) break;
            ng[group_of_node[u]]=0;
        }
    }
    return ans;
}

std::vector<int> intersect(const std::vector<std::array<int,2> > gu, const std::vector<int> ci, const int minv) {
    auto i = gu.begin();
    auto ie = gu.end();
    auto j = ci.begin();
    auto je = ci.end();
    while(i!=ie && i->at(0)<=minv) ++i;
    while(j!=je && *j<=minv) ++j;
    std::vector<int> ans;
    while(i!=ie && j!=je) {
        if(i->at(0)<*j) ++i;
        else if(*j<i->at(0)) ++j;
        else {
            ans.push_back(*j);
            ++i;++j;
        }
    }
    return ans;
}

// TODO: separate ccs
// TODO: get a primal bound even when solution is not integer
std::vector<mcolumn> columns; // TODO: fill this with columns and pass them to children
int primal_bound=0;
void bnp(mgraph mg) {
    Master m(mg, columns_for(columns, mg));
    m.setup_highs();
    // m.hot_start();
    double dv = m.run();

    int dual_bound = int(0.001+dv);
    if(dual_bound<=primal_bound) return;
    const std::vector<double>& cv = m.highs.getSolution().col_value;
    std::cerr<<cv.size()<<' '<<m.columns.size()<<std::endl;
    assert(cv.size()==m.columns.size());

    // find where to split
    double mincv=1;
    std::map<std::array<int,2>, double> tmp;
    for(int i=0; i<(int)cv.size(); ++i) if(cv[i]>0.0001) {
        if(cv[i]<mincv) mincv=cv[i];
        std::cerr<<"col "<<i<<" "<<cv[i]<<"*";
        std::cerr<<m.columns[i].value<<'\t';
        for(auto &u: m.columns[i].nodes) {
            std::cerr<<u<<' ';
            for(auto &v: intersect(m.g.g[u], m.columns[i].nodes, u+1)) {
                tmp[{u,v}]+=cv[i];
            }
        }
        std::cerr<<'\n';
    }
    double md=2;
    std::array<int,2> split={0,0};
    for(auto &[uv,w]: tmp) {
        if(fabs(w-0.5)<md) {
            md=fabs(w-0.5);
            split=uv;
        }
    }
    if(mincv>0.999) { // solution is integer, found a primal bound
        std::cerr<<"BNP integer: "<<dv<<std::endl;
        primal_bound = dual_bound;
        // TODO: save best solution //best_solution=mg;
    } else { // split
        std::cerr<<"BNP frac: "<<dv<<" split on:"<<split[0]<<' '<<split[1]<<std::endl;
        mgraph mg1(mg);
        mgraph mg2(mg);
        mg1.merge_edge(split[0],split[1]);
        mg2.remove_edge(split[0], split[1]);
        bnp(mg1);
        bnp(mg2);
    }
}

void bnp_main(problem p) {
    mgraph mg0(p);
	mg0.simplify();
    bnp(mg0);
}
