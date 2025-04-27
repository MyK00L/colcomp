#pragma once

#include <vector>
#include <array>
#include "mgraph.cpp"
#include <queue>

struct pricing_problem {
    mgraph g;
    std::vector<double> cost;
    bool is_node_useless(int u, const std::vector<bool>& nu, const std::vector<std::vector<bool> >& eu) {
        std::vector<int> ma(g.C(),0);
        for(auto &[v,w]:g.g[u]) if(!nu[v] && !eu[u][v]) ma[g.c[v][0]]=std::max(ma[g.c[v][0]],w);
        double x = -cost[u];
        for(auto &y:ma)x+=y;
        return x<=0;
    }
    bool is_edge_useless(int u, int v, const std::vector<bool>& nu, const std::vector<std::vector<bool> >& eu) {
        std::vector<int> ma(g.C(),0);
        for(auto &[j,w]:g.g[u]) if(!nu[j] && ! eu[u][j]) ma[g.c[j][0]]=std::max(ma[g.c[j][0]],w);
        for(auto &[j,w]:g.g[v]) if(!nu[j] && ! eu[u][j]) ma[g.c[j][0]]=std::max(ma[g.c[j][0]],w);
        int uvw=0;
        for(auto &[j,w]: g.g[u]) if(j==v) uvw=w;
        for(auto &cu: g.c[u]) ma[cu]=0;
        for(auto &cv: g.c[v]) ma[cv]=0;
        double x = -cost[u]-cost[v]+uvw;
        for(auto &y:ma) x+=y;
        return x<=0;
    }
    void preprocess() {
        std::vector<bool> nu(g.n(),0); // 1 if node is useless
        std::vector<std::vector<bool> > eu(g.n(), std::vector<bool>(g.n(),0)); // 1 if edge is useless
        std::vector<bool> inq(g.n(),1);
        std::queue<int> q;
        for(int i=0; i<g.n(); ++i) q.push(i);
        int nuseless=0;
        bool upd=1;
        while(upd) {
            while(!q.empty()) {
                int u = q.front();
                q.pop();
                if(is_node_useless(u,nu,eu)) {
                    nu[u]=1;
                    ++nuseless;
                    std::vector<std::array<int,2> > tredg;
                    for(auto &[v,_w]: g.g[u]) {
                        eu[u][v] = eu[v][u] = 1;
                        tredg.push_back({u,v});
                        if(!nu[v] && !inq[v]) {
                            inq[v]=1;
                            q.push(v);
                        }
                    }
                }
                inq[u]=0;
            }
            upd=0;
            for(int u=0; u<g.n() && !upd; ++u) if(!nu[u]) for(auto &[v,w]:g.g[u]) if(!nu[v] && !eu[u][v]) {
                if(is_edge_useless(u, v, nu, eu)) {
                    upd=1;
                    eu[u][v]=eu[v][u]=1;
                    inq[u]=1;
                    inq[v]=1;
                    q.push(u);
                    q.push(v);
                    break;
                }
            }
        }
        for(int i=0; i<g.n(); ++i) {
            for(auto &[j,w]: g.g[i]) if(eu[i][j]) w=-g.orig_n*g.orig_n-1;
            if(nu[i]) cost[i] += g.orig_n+1;
        }

        int e0=0,e1=0;
        for(int i=0; i<g.n(); ++i) for(auto &[j,w]: g.g[i]) if(j>i) {
            e0++;
            e1+=!eu[i][j];
        }
        // std::cerr<<"|V| |E|: "<<g.n()<<' '<<e0<<"  ->  "<<g.n()-nuseless<<' '<<e1<<std::endl;
    }
    pricing_problem(mgraph _g, std::vector<double> duals): g(_g), cost(duals) {
        preprocess();
    }
};
