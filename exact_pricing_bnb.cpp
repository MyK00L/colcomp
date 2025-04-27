#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <queue>
#include <array>
#include "common.cpp"
#include "pricing.hpp"

struct epbnb_state {
    bitvec nodes; // nodes taken
    bitvec untakeable; // nodes discarded or taken or with uncompatible colors
    double score; // score achieved until now
    epbnb_state(const pricing_problem& p): nodes(p.g.n()), untakeable(p.g.n()), score(0) {}
    void take(const int u, const pricing_problem& p) {
        score-=p.cost[u];
        for(auto &[v,w]: p.g.g[u]) if(nodes[v]) score+=w;
        nodes.set1(u);
        untakeable.set1(u);
        for(int v=0; v<p.g.n(); ++v) if(!untakeable[v]) {
            if(!compatible(p.g.c[u],p.g.c[v])) untakeable.set1(v);
            // for(auto &cv: p.g.c[v]) if(colors[cv]) untakeable.set1(v);
        }
        // remove dummy node from taken
        nodes.set0(p.g.n()-1);
    }
    void forbid(const int u) {
        untakeable.set1(u);
    }
    std::array<epbnb_state,2> branch(const pricing_problem& p) const {
        int minv=p.g.n();
        for(int u=0; u<p.g.n()&&minv==p.g.n(); ++u) if(nodes[u]) {
            for(auto &[v,w]:p.g.g[u]) if(!untakeable[v]) {
                minv = v;
                break;
            }
        }
        epbnb_state take_node(*this);
        epbnb_state forbid_node(*this);
        take_node.take(minv, p);
        forbid_node.forbid(minv);
        return {take_node,forbid_node};
    }
    double upper_bound(const pricing_problem& p) const {

        int cl = 0; // colors left to take
        {
            bitvec cvis(p.g.C());
            for(int i=0; i<p.g.n(); ++i) if(!untakeable[i] && !cvis[p.g.c[i][0]]) {
                cl+=1;
                cvis.set1(p.g.c[i][0]);
            }
        }
        bitvec reach(p.g.n());
        while(cl) {
            int newcl=0;
            bitvec newreach(p.g.n());
            bitvec cvis(p.g.C());
            std::queue<std::array<int,2> > q;
            for(int u=0; u<p.g.n(); ++u) if(nodes[u]) {
                for(auto &[v,w]:p.g.g[u]) if(!untakeable[v] && !newreach[v]) {
                    q.push({v,1});
                    newreach.set1(v);
                    if(!cvis[p.g.c[v][0]]) {
                        cvis.set1(p.g.c[v][0]);
                        newcl+=1;
                    }
                }
            }
            while(!q.empty()) {
                auto [u,d] = q.front();
                q.pop();
                if(d==cl) continue;
                for(auto &[v,w]:p.g.g[u]) if(!untakeable[v] && !newreach[v]) {
                    q.push({v,d+1});
                    newreach.set1(v);
                    if(!cvis[p.g.c[v][0]]) {
                        cvis.set1(p.g.c[v][0]);
                        newcl+=1;
                    }
                }
            }
            reach=newreach;
            assert(newcl<=cl);
            if(newcl==cl) break;
            cl=newcl;
        }
        if(cl==0) return score;
        std::vector<double> bbc(p.g.C(),0); // best reachable value by color
        std::vector<int> bc(p.g.C(), 0); // best connection by color, for internal loop
        for(int u=0; u<p.g.n(); ++u) if(reach[u]) {
            fill(bc.begin(),bc.end(),0);
            for(auto &[v,w]: p.g.g[u]) {
                if(nodes[v]) {
                    bc[p.g.c[v][0]] = w;
                } else if (reach[v]) {
                    bc[p.g.c[v][0]] = std::max(bc[p.g.c[v][0]],w); // TODO: better bound
                }
            }
            double us=-p.cost[u];
            for(auto &sc:bc) us+=sc;
            if(us>bbc[p.g.c[u][0]]) bbc[p.g.c[u][0]]=us;
        }
        double ans = score;
        for(auto &sc: bbc) ans+=sc;
        return ans;
    }
};

void epbnb(const pricing_problem& p, const epbnb_state state, const double dual, double& primal, bitvec& ans) {
    if(state.score>primal) {
        primal=state.score;
        ans=state.nodes;
    }
    if(primal>=dual) return;
    auto [s1,s2] = state.branch(p);
    const double d1 = std::min(dual,s1.upper_bound(p));
    const double d2 = std::min(dual,s2.upper_bound(p));
    if(d2>d1) {
        epbnb(p,s2,d2,primal,ans);
        epbnb(p,s1,d1,primal,ans);
    } else {
        epbnb(p,s1,d1,primal,ans);
        epbnb(p,s2,d2,primal,ans);
    }
}

std::vector<mcolumn> exact_pricing_bnb(pricing_problem p) {
    // add dummy node
    p.g.g.push_back({});
    p.g.c.push_back({});
    p.g.groups.push_back({});
    p.cost.push_back(0);
    for(int v=0; v<p.g.n()-1; ++v) {
        p.g.g.back().push_back({v,0});
    }
    epbnb_state s0(p);
    s0.nodes.set1(p.g.n()-1);
    s0.untakeable.set1(p.g.n()-1);
    double primal=0;
    double dual = p.g.n()*p.g.n();
    bitvec ans = s0.nodes;
    epbnb(p,s0,dual,primal,ans);
    if(primal>1e-6) {
        mcolumn c;
        double val = primal;
        for(int i=0; i<p.g.n()-1; ++i) if(ans[i]) {
            c.nodes.push_back(i);
            val+=p.cost[i];
        }
        c.value = int(round(val));
        return {c};
    }
    return {};
}

std::vector<mcolumn> hot_start_cols(mgraph mg) {
    pricing_problem p(mg,std::vector<double>(mg.n(),0));
    // add dummy node
    p.g.g.push_back({});
    p.g.c.push_back({});
    p.g.groups.push_back({});
    p.cost.push_back(0);
    for(int v=0; v<p.g.n()-1; ++v) {
        p.g.g.back().push_back({v,0});
    }
    epbnb_state s0(p);
    s0.nodes.set1(p.g.n()-1);
    s0.untakeable.set1(p.g.n()-1);

    for(int i=0; i<p.g.n(); ++i) assert(p.cost[i]==0);
    std::vector<mcolumn> res;
    while(1) {
        epbnb_state si = s0;
        double primal=0;
        double dual = si.upper_bound(p);
        bitvec ans(p.g.n());
        epbnb(p,si,dual,primal,ans);

        if(primal<0.5) break;
        mcolumn c;
        for(int u=0; u<p.g.n()-1; ++u) if(ans[u]) {
            c.nodes.push_back(u);
            s0.forbid(u);
        }
        c.value = int(round(primal));
        res.push_back(c);
    } // TODO: check duplicates
    return res;
}
