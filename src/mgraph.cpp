# pragma once
#include <cmath>
#include <iostream>
#include <iterator>
#include <queue>
#include <vector>
#include <array>
#include <algorithm>
#include "common.cpp"

/// graph where nodes can be merged and have colors
struct mgraph {
    /// original number of nodes in the whole graph
    int orig_n;
    /// number of colors in this component
	int _C;
    /// original nodes grouped in each node
    std::vector<std::vector<int> > groups;
    /// colors of each node (may be remapped)
    std::vector<std::vector<int> > c;
    /// edges of each node
    std::vector<std::vector<std::array<int,2> > > g;
    mgraph(const problem& p) : orig_n((int)p.n()), _C(0), groups(p.n()), c(p.n()), g(p.n()) {
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
	mgraph(): orig_n(0), _C(0) {}
    /// returns the size of the graph
    int n() const {
        return (int)c.size();
    }
    /// returns the number of colors
    int C() const {
        return _C;
    }
	/// returns the number of edges in O(m)
	int m() const {
		int ans = 0;
		for(auto &gu:g) for(auto &[_v,w]:gu) ans+=w;
		return ans/2;
	}
	/// returns indeces original nodes that are in this component
	std::vector<int> orig_nodes() const {
		std::vector<int> ans;
		for(auto &grp: groups) for(auto &ou: grp) ans.push_back(ou);
		return ans;
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
	/// remaps colors in [0,#colors)
	void simplify_colors() {
		// TODO
	}
	/// applies rule 2 to remove edges of some cuts
	void rule2() {
		// TODO
	}
	bool is_colorful() const {
		std::vector<bool> visc(C(),0);
		for(auto &cu: c) {
			for(auto &cui: cu) {
				if(visc[cui]) return 0;
				visc[cui]=1;
			}
		}
		return 1;
	}
    /// returns true if all connected components are colorful
    bool is_trivial() const {
        std::vector<bool> visn(n(), 0);
        for(int i=0; i<n(); ++i) if(!visn[i]) {
            visn[i] = 1;
            std::vector<bool> visc(C(), 0);
            std::queue<int> q;
            q.push(i);
            while(!q.empty()) {
                int u = q.front();
                q.pop();
                for(auto &cu: c[u]) {
                    if(visc[cu]) return 0;
                    visc[cu] = 1;
                }
                for(auto &[v,_w]: g[u]) if(!visn[v]) {
                    visn[v]=1;
                    q.push(v);
                }
            }
        }
        return 1;
    }
    /// splits the graph into connected components
    std::vector<mgraph> split_components() const {
		std::vector<std::vector<int> > ccs;
        std::vector<bool> vis(n(),0);
        for(int u0=0; u0<n(); ++u0) if(!vis[u0]) {
            ccs.push_back({});
			std::queue<int> q;
            q.push(u0);
            vis[u0]=1;
            while(!q.empty()) {
                int u = q.front();
                q.pop();
				ccs.back().push_back(u);
                for(auto &[v,_w]: g[u]) if(!vis[v]) {
                    vis[v]=1;
                    q.push(v);
                }
            }
			sort(ccs.back().begin(),ccs.back().end());
    	}
		std::vector<int> newi(n());
		std::vector<mgraph> ans;
		for(auto& cc:ccs) {
			mgraph mg;
			int newn=0;
			for(auto &u:cc) newi[u]=newn++;
			mg.g.resize(newn);
			mg.c.resize(newn);
			mg.groups.resize(newn);
			mg.orig_n=orig_n;
			mg._C = _C;
			for(auto &u:cc) {
				mg.c[newi[u]] = c[u];
				mg.g[newi[u]].reserve(g[u].size());
				for(auto &[v,w]:g[u]) {
					mg.g[newi[u]].push_back({newi[v],w});
				}
				mg.groups[newi[u]] = groups[u];
			}
			mg.simplify_colors();
			ans.push_back(mg);
		}
		return ans;
	}
    /// merges edge between u and v into a new node at the end of the graph
    /// returns the number of edges between u and v
    int merge_edge(int u, int v) {
        if(u>v) std::swap(u,v);
        const int ret = edge_val(u, v);
        
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

        // simplify();

        return ret;
    }
    /// remove edge (u,v) and returns its weight
    int remove_edge(int u, int v) {
        int ret = edge_val(u, v);
        auto ituv = std::lower_bound(g[u].begin(),g[u].end(),std::array<int,2>({v,0}),[](const auto&a, const auto&b){return a[0]<b[0];});
        if(ituv!=g[u].end() && ituv->at(0)==v) g[u].erase(ituv);
        auto itvu = std::lower_bound(g[v].begin(),g[v].end(),std::array<int,2>({u,0}),[](const auto&a, const auto&b){return a[0]<b[0];});
        if(itvu!=g[v].end() && itvu->at(0)==u) g[v].erase(itvu);
        // simplify();
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
    int tight_partition_bound() const {
        std::vector<int> nc;
        for(auto&cs: c) for(auto &ci: cs) {
            if(ci >= (int)nc.size()) nc.resize(ci+1,0);
            nc[ci]+=1;
        }
        int max_edges = 0;
        for(auto &c: nc) {
            max_edges += c*(c-1)/2;
        }
        return std::min((int)n(), int(sqrt(max_edges*2+0.25)+0.5+0.001));
    }
    void debug_print() const {
        std::cerr<<"mgraph:"<<std::endl;
        for(int u=0; u<n(); ++u) {
            std::cerr<<u<<"(";
            for(auto &uc: c[u])std::cerr<<uc<<",";
            std::cerr<<")[";
            for(auto &uo: groups[u])std::cerr<<uo<<",";
            std::cerr<<"]\n";
        }
        for(int u=0; u<n(); ++u) for(const auto &[v,w]: g[u]) {
            if(v>u) {
                std::cerr<<u<<"-"<<v<<" "<<w<<std::endl;
            }
        }
    }
};

