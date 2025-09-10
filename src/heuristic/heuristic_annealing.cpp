#pragma once

#include <cassert>
#include <cmath>
#include <iterator>
#include <set>
#include <vector>
#include <algorithm>
#include <iostream>
#include "common.cpp"
#include "heuristic_beamsearch.cpp"

struct annealing_node {
    struct undo_info {
        int i;
        int ii;
        int j;
        int jj;
    };
    int score;
    static std::vector<int> color;
    static std::set<std::array<int,2> > edges;
    std::vector<std::vector<int> > groups; // groups[i]= nodes in group i
    /// returns number of edges between u and elements of groups[i]
    int calc(int u, int i) {
        int ans=0;
        for(int j=0; j<(int)groups[i].size(); ++j)
            ans+=edges.count({std::min(u,groups[i][j]),std::max(u,groups[i][j])});
        return ans;
    }
    /// moves node groups[i][ii] to group j
    undo_info movec(int i, int ii, int j) {
        if(j==(int)groups.size()) {
            groups.push_back({});
        }
        if(groups[i].size()==1) {
            groups[i].swap(groups.back());
            if(j==(int)groups.size()-1) j=i;
            i=groups.size()-1;
        }
        std::swap(groups[i][ii],groups[i].back());
        int x = groups[i].back();
        groups[i].pop_back();
        score-=calc(x, i);
        score+=calc(x, j);
        groups[j].push_back(x);
        if(groups.back().size()==0) groups.pop_back();
        return undo_info{i,ii,j,-1};
    }
    /// swap nodes group[i][ii] and group[j][jj]
    undo_info swapc(int i, int ii, int j, int jj) {
        score-=calc(groups[i][ii],i);
        score-=calc(groups[j][jj],j);
        std::swap(groups[i][ii], groups[j][jj]);
        score+=calc(groups[i][ii],i);
        score+=calc(groups[j][jj],j);
        return undo_info{i,ii,j,jj};
    }
    void undo(const undo_info& info) {
        const auto [i,ii,j,jj]=info;
        if(jj==-1) { // move
            movec(j,groups[j].size()-1,i);
        } else { // swap
            swapc(i,ii,j,jj);
        }
    }
    // neighbourhood:
    // - swap 2 nodes of the same color
    // - move a node to another group that does not have its color
    /// samples a neighbourhood, moving to a random neighbour, returns info to undo
    undo_info sample_neighbourhood(Rng& rng) {
        bool sn = 0;
        int i = rng()%groups.size();
        int j = rng()%groups.size(); if(j==i) ++j;
        if(j==(int)groups.size()) groups.push_back({});

        std::vector<int> ci,cj;
        for(auto &e:groups[i]) ci.push_back(color[e]);
        for(auto &e:groups[j]) cj.push_back(color[e]);
        std::sort(ci.begin(),ci.end());
        std::sort(cj.begin(),cj.end());
        std::vector<int> cij;
        std::set_intersection(ci.begin(),ci.end(),cj.begin(),cj.end(),std::back_inserter(cij));
        if(ci==cj) sn=1;
        else if(cij.size()==0) sn=0;
        else sn=rng()%2;
        if(sn) { // swap neighbourhood
            int c = cij[rng()%cij.size()];
            int ii=0,jj=0;
            while(color[groups[i][ii]]!=c) ++ii;
            while(color[groups[j][jj]]!=c) ++jj;
            return swapc(i,ii,j,jj);
        } else { // move neighbourhood
            std::vector<int> un;
            std::set_union(ci.begin(),ci.end(),cj.begin(),cj.end(),std::back_inserter(un));
            std::vector<int> df;
            std::set_difference(un.begin(),un.end(),cij.begin(),cij.end(),std::back_inserter(df));
            int c = df[rng()%df.size()];
            if(!std::binary_search(ci.begin(),ci.end(),c)) std::swap(i,j);
            int ii=0;
            while(color[groups[i][ii]]!=c) ++ii;
            return movec(i,ii,j);
        }
    }
};
std::vector<int> annealing_node::color = {};
std::set<std::array<int,2> > annealing_node::edges = {};

void simulated_annealing(problem p) {
    annealing_node::color.clear();
    annealing_node::edges.clear();
    for(auto &e:p.edges) {
        int u = std::min(e[0],e[1]);
        int v = std::max(e[0],e[1]);
        annealing_node::edges.insert({u,v});
    }
    annealing_node::color.reserve(p.color.size());
    for(auto &c: p.color) annealing_node::color.push_back(c);
    std::vector<std::vector<int> > g0 = beamsearch(p,1);
    annealing_node node{0,g0};
    annealing_node best = node;
    double temp=10;//p.edges.size();
    Rng rng(4);
    int it=0;
    while(temp>0.0001) {
        int old_score = node.score;
        auto ui = node.sample_neighbourhood(rng);
        if(node.score>best.score) best = node;
        if(it++%100==0)std::cerr<<temp<<' '<<old_score<<" -> "<<node.score<<std::endl;
        if(node.score<old_score && std::exp((node.score-old_score)/temp)<rng.f01()) {
            node.undo(ui);
            assert(node.score==old_score);
        }
        temp*=0.99;
    }
    std::cerr<<"annealing: "<<best.score<<std::endl;
}
