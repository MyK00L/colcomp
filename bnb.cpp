#include "common.cpp"
#include "mgraph.cpp"
#include <algorithm>
#include <stack>

struct bnbstate {
    bnbstate(mgraph _mg): mg(_mg), dual_bound(0) {
        std::vector<int> cfreq(1,0);
        dual_bound=mg.value;
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
    std::array<bnbstate,2> split() const {
        mgraph mg1(mg),mg2(mg);
        int u=0; // TODO: better choice for split
        int v=mg.g[u][0][0];
        mg1.merge_edge(u,v);
        mg2.remove_edge(u,v);
        return {bnbstate(mg1), bnbstate(mg2)};
    }
    bool over() const {
        return mg.n()==0;
    }
};

void bnb(problem p) {
    int primal_bound = 0;
    mgraph best(p);
    std::stack<bnbstate> st;
    st.push(bnbstate(mgraph(p)));
    while(!st.empty()) {
        auto top = st.top();
        st.pop();
        auto t = top.split();
        if(t[0].dual_bound<t[1].dual_bound) std::swap(t[0],t[1]);
        for(auto &o: t) if(o.dual_bound>primal_bound) {
            if(o.over()) {
                primal_bound = o.dual_bound;
                best = o.mg;
            } else {
                st.push(o);
            }
        }
    }
    std::cerr<<"bnb: "<<primal_bound<<std::endl;
}
