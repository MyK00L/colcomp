#pragma once

#include <cmath>
#include <iterator>
#include <limits>
#include <numeric>
#include <vector>
#include <array>
#include <algorithm>
#include "pricing.hpp"
#include "pricing.hpp"

template<class T>
struct maxheap {
    std::vector<size_t> pos;
    std::vector<std::pair<int,T> > heap;
	maxheap(size_t n): pos(n), heap(n+1) {
        std::iota(pos.begin(),pos.end(),1);
        for(size_t i=1; i<=n; ++i) heap[i]={int(i)-1,T()};
    }
    void update(int _i, T val) {
        size_t i = pos[_i];
        heap[i].second=val;
        // push up
        while(i>1 && heap[i].second>heap[i/2].second) {
            swap(heap[i],heap[i/2]);
            pos[heap[i].first]=int(i);
            pos[heap[i/2].first]=int(i/2);
            i=i/2;
        }
        // push down
        while(i*2<heap.size()) {
            size_t j = i*2+1;
            if(i*2+1>=heap.size() || heap[i*2].second>heap[i*2+1].second) j=i*2;
            if(heap[j].second>heap[i].second) {
                swap(heap[i],heap[j]);
                pos[heap[i].first]=i;
                pos[heap[j].first]=j;
                i=j;
            } else {
                break;
            }
        }
    }
    std::pair<int,T> top() {
        return heap[1];
    }
};

void filter_edges(std::vector<std::array<int,2> >& edges, int min) {
    edges.erase(std::remove_if(edges.begin(),edges.end(),[&min](std::array<int,2> e){return e[0]<min;}),edges.end());
}
void sort_by_value(std::vector<std::array<int,2> >& edges, const std::vector<double>& duals){
    std::sort(edges.begin(),edges.end(),[&duals](std::array<int,2>& a, std::array<int,2>& b){return a[1]-duals[a[0]]<b[1]-duals[b[0]];});
}
struct subgraph {
    double value;
    std::vector<int> nodes;
    std::vector<int> uc; // used colors
    std::vector<std::array<int,2> > edges;
    int min_node() const {
        return nodes[0];
    }
    void upd_nxt(const pricing_problem& pp) {
        while(edges.size() && !compatible(pp.g.c[edges.back()[0]], uc)) edges.pop_back();
    }
    double nxt_val(const std::vector<double>& duals) const {
        if(edges.empty()) return std::numeric_limits<double>::min();
        return value+edges.back()[1]-duals[edges.back()[0]];
    }
    size_t nxt_ind() const {
        return edges.back()[0];
    }
};
subgraph single_node_subgraph(const pricing_problem& pp, const std::vector<double>& duals, int i) {
    double value = -duals[i];
    std::vector<int> nodes = {i};
    std::vector<int> uc = pp.g.c[i];
    std::vector<std::array<int,2> > edges = pp.g.g[i];
    return subgraph{value,nodes,uc,edges};
}
/// assumes b is a subgraph for a single node
/// does not remove untakable edges and does not update value
subgraph merge(subgraph& a, const subgraph& b, const std::vector<double>& duals) {
    std::vector<int> nodes;
    nodes.reserve(a.nodes.size()+b.nodes.size());
    std::merge(a.nodes.begin(), a.nodes.end(), b.nodes.begin(), b.nodes.end(), std::back_inserter(nodes));
   
    std::vector<int> uc;
    uc.reserve(a.uc.size()+b.uc.size());
    std::merge(a.uc.begin(), a.uc.end(), b.uc.begin(), b.uc.end(), std::back_inserter(uc));

    std::vector<std::array<int,2> > edges;
    std::sort(a.edges.begin(),a.edges.end());
    edges.reserve(a.edges.size()+b.edges.size());
    size_t ai=0,bi=0;
    while(ai<a.edges.size() && a.edges[ai][0]<a.min_node()) ++ai;
    while(bi<b.edges.size() && b.edges[bi][0]<a.min_node()) ++bi;
    while(ai<a.edges.size() && bi<b.edges.size()) {
        if(a.edges[ai][0]<b.edges[bi][0]) {
            edges.push_back(a.edges[ai]);
            ++ai;
        } else if (a.edges[ai][0]>b.edges[bi][0]) {
            edges.push_back(b.edges[bi]);
            ++bi;
        } else {
            edges.push_back({a.edges[ai][0],a.edges[ai][1]+b.edges[bi][1]});
            ai++;
            bi++;
        }
    }
    while(ai<a.edges.size()) {
        edges.push_back(a.edges[ai]);
        ++ai;
    }
    while(bi<b.edges.size()) {
        edges.push_back(b.edges[bi]);
        ++bi;
    }

    sort_by_value(a.edges, duals);
    sort_by_value(edges, duals);

    return subgraph{0.0,nodes,uc,edges};
}

// keepa all edges and check for colors incrementally
// only choose nodes bigger than the first one to avoid doubles
/// beam search O(extneigh * beam_width * depth)
std::vector<column> heuristic_pricing_beamsearch(const pricing_problem& pp) {
    const std::vector<double>& duals = pp.cost;
    const size_t beam_width = 200;
    std::vector<subgraph> singles;
    singles.reserve(pp.g.n());
    for(int i=0; i<int(pp.g.n()); ++i) {
        singles.push_back(single_node_subgraph(pp, duals, i));
    }
    std::vector<subgraph> beam = singles;
    double max_v = 0;
    std::vector<int> max_nz;
    do {
        //std::cerr<<beam.size()<<' '<<m.n()<<std::endl;
        maxheap<double> heap(beam.size());
        std::vector<subgraph> next_beam;
        next_beam.reserve(beam_width);
        for(int i=0; i<(int)beam.size(); ++i) heap.update(i, beam[i].nxt_val(duals));
        while(next_beam.size()<beam_width && heap.top().second!=std::numeric_limits<double>::min()) {
            auto [i,v] = heap.top();
            subgraph nb = merge(beam[i],singles[beam[i].nxt_ind()], duals);
            nb.value = v;
            if(v>max_v) {
                max_v=v;
                max_nz = nb.nodes;
            }
            beam[i].edges.pop_back();
            beam[i].upd_nxt(pp);
            heap.update(i,beam[i].nxt_val(duals));
            nb.upd_nxt(pp);
            next_beam.push_back(nb);
        }
        beam=next_beam;
    } while(!beam.empty());

    //std::cerr<<"max_v: "<<max_v<<" generated column: ";
    //for(size_t i=0; i<max_nz.size(); ++i) std::cerr<<max_nz[i]<<' ';
    //std::cerr<<std::endl;
    //for(size_t i=0; i<max_nz.size(); ++i) std::cerr<<duals[max_nz[i]]<<' ';
    //std::cerr<<std::endl;
    if(max_v>1e-3) return {column{max_v,max_nz}};
    return {};
}
