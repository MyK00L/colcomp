# pragma once
#include <iterator>
#include <queue>
#include <vector>
#include <array>
#include <algorithm>
#include "common.cpp"

#include <iostream>

/// graph where nodes can be merged and have colors
struct mgraph {
    /// original n
    int orig_n;
    /// number of edges taken until now
    int value;
    int _C;
    std::vector<std::vector<int> > final_groups;
    /// original nodes grouped in each node
    std::vector<std::vector<int> > groups;
    /// colors of each node
    std::vector<std::vector<int> > c;
    /// edges of each node
    std::vector<std::vector<std::array<int,2> > > g;
    mgraph(problem p) : orig_n((int)p.n()), value(0), _C(0), final_groups(0), groups(p.n()), c(p.n()), g(p.n()) {
        for(auto [u,v] : p.edges) if(p.color[u]!=p.color[v]) {
            g[u].push_back({(int)v,1});
            g[v].push_back({(int)u,1});
        }
        for(int i=0; i<(int)g.size(); ++i) groups[i].push_back(i);
        for(size_t i=0; i<g.size(); ++i) sort(g[i].begin(),g[i].end());
        for(size_t i=0; i<p.color.size(); ++i) {
            c[i].push_back(p.color[i]);
            _C=std::max(_C, (int)p.color[i]+1);
        }
        assert(g.size()==c.size());
    }
    /// returns the size of the graph
    int n() const {
        return (int)c.size();
    }
    /// returns the number of colors
    int C() const {
        return _C;
    }
    /// returns the weight of the edge between u and v, 0 if there is no such edge
    int edge_val(int u, int v) const {
        auto it = std::lower_bound(g[u].begin(),g[u].end(),std::array<int,2>({v,0}),[](const auto&a, const auto&b){return a[0]<b[0];});
        if(it==g[u].end() || it->at(0)!=v) return 0;
        return it->at(1);
    }
    /// returns true if u and v have no colors in common and can thus be merged
    bool can_merge(int u, int v) const {
        return compatible(c[u], c[v]);
    }
    void delete_nodes(std::vector<int> v) {
        sort(v.begin(),v.end());
        int N = n();
        auto vi = v.begin();
        int ni=0;
        std::vector<int> newi(N,-1);
        for(int i=0; i<N; ++i) {
            if(vi!=v.end() && *vi==i) {++vi; continue;}
            newi[i]=ni;
            if(ni<i) {
                g[ni]=std::move(g[i]);
                c[ni]=std::move(c[i]);
                groups[ni]=std::move(groups[i]);
            }
            ++ni;
        }
        g.erase(g.begin()+ni,g.end());
        c.erase(c.begin()+ni, c.end());
        groups.erase(groups.begin()+ni, groups.end());
        for(auto &gi: g) {
            N=(int)gi.size();
            ni=0;
            for(int i=0; i<N; ++i) {
                if(newi[gi[i][0]]==-1) continue;
                gi[i][0]=newi[gi[i][0]];
                if(ni<i) {
                    gi[ni]=std::move(gi[i]);
                }
                ++ni;
            }
            gi.erase(gi.begin()+ni,gi.end());
        }
    }
    /// removes components that can be completely merged
    void simplify() {
        std::vector<bool> vis(n(),0);
        std::vector<int> to_delete;
        for(int u0=0; u0<n(); ++u0) if(!vis[u0]) {
            std::queue<int> q;
            q.push(u0);
            vis[u0]=1;
            std::vector<bool> visc(n(),0);
            std::vector<int> nodes;
            bool del=1;
            while(!q.empty()) {
                int u = q.front();
                q.pop();
                nodes.push_back(u);
                for(auto &[v,_w]: g[u]) if(!vis[v]) {
                    vis[v]=1;
                    q.push(v);
                }
                for(auto col: c[u]) {
                    if(col>=(int)visc.size()) visc.resize(col+1,0);
                    if(visc[col]) del=0;
                    visc[col]=1;
                }
            }
            if(del==1) {
                std::vector<int> ng;
                for(auto u: nodes) {
                    to_delete.push_back(u);
                    for(auto &[v,w]: g[u]) if(v>u) value+=w;
                    std::vector<int> tng;
                    tng.swap(ng);
                    std::merge(std::move_iterator(tng.begin()),std::move_iterator(tng.end()),std::move_iterator(groups[u].begin()),std::move_iterator(groups[u].end()),std::back_insert_iterator(ng));
                }
                final_groups.push_back(ng);
            }
        }
        delete_nodes(to_delete);
    }
    /// merges edge between u and v into a new node at the end of the graph
    /// returns the number of edges between u and v
    int merge_edge(int u, int v) {
        if(u>v) std::swap(u,v);
        const int ret = edge_val(u, v);
        value+=ret;
        std::sort(g[u].begin(),g[u].end());
        std::sort(g[v].begin(),g[v].end());
        auto pu = std::move_iterator(g[u].begin());
        const auto eu = std::move_iterator(g[u].end());
        auto pv = std::move_iterator(g[v].begin());
        const auto ev = std::move_iterator(g[v].end());
        std::vector<std::array<int,2> > ngi; ngi.reserve(g[u].size()+g[v].size());
        while(pu!=eu || pv!=ev) {
            if(pu==eu || (pv!=ev && pv->at(0)<pu->at(0))) {
                if(compatible(c[u],c[pv->at(0)])) ngi.push_back(*pv);
                ++pv;
            } else if(pv==ev || (pu!=eu && pu->at(0)<pv->at(0))) {
                if(compatible(c[v],c[pu->at(0)])) ngi.push_back(*pu);
                ++pu;
            } else {
                ngi.push_back({pu->at(0),pu->at(1)+pv->at(1)});
                ++pu;
                ++pv;
            }
        }
        g.push_back(std::move(ngi));
        for(auto [j,w]: g.back()) g[j].push_back({(int)g.size()-1,w});

        std::vector<int> nc; nc.reserve(c[u].size()+c[v].size());
        std::merge(std::move_iterator(c[u].begin()),std::move_iterator(c[u].end()),std::move_iterator(c[v].begin()),std::move_iterator(c[v].end()),std::back_insert_iterator(nc));
        c.push_back(std::move(nc));

        std::vector<int> ng; ng.reserve(groups[u].size()+groups[v].size());
        std::merge(std::move_iterator(groups[u].begin()), std::move_iterator(groups[u].end()), std::move_iterator(groups[v].begin()), std::move_iterator(groups[v].end()), std::back_insert_iterator(ng));
        groups.push_back(std::move(ng));

        delete_nodes({u,v});

        simplify();

        return ret;
    }
    /// remove edge (u,v) and returns its weight
    int remove_edge(int u, int v) {
        int ret = edge_val(u, v);
        auto ituv = std::lower_bound(g[u].begin(),g[u].end(),std::array<int,2>({v,0}),[](const auto&a, const auto&b){return a[0]<b[0];});
        if(ituv!=g[u].end() && ituv->at(0)==v) g[u].erase(ituv);
        auto itvu = std::lower_bound(g[v].begin(),g[v].end(),std::array<int,2>({u,0}),[](const auto&a, const auto&b){return a[0]<b[0];});
        if(itvu!=g[v].end() && itvu->at(0)==u) g[v].erase(itvu);
        simplify();
        return ret;
    }
    bool has_edge(int u, int v) const {
        if(g[u].size()>g[v].size()) return has_edge(v, u);
        return std::binary_search(g[u].begin(),g[u].end(),std::array<int,2>{v,0},[](const std::array<int,2>& a, const std::array<int,2>& b){
            return a[0]<b[0];
        });
    }
    int value_for_col(const std::vector<int>& col) const {
        int wsum=0;
        for(int i=0; i<(int)col.size(); ++i) {
            for(int j=i+1; j<(int)col.size(); ++j) {
                wsum += edge_val(col[i],col[j]);
            }
        }

        return wsum;
    }
};

