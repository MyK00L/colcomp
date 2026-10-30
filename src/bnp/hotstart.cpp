#include <iterator>
#include <set>

#include "common.cpp"
#include "heuristic/heuristic_karger.cpp"

void gen_with_greedy(const problem &p, const std::vector<std::vector<int>> &g, std::vector<std::vector<int>>& cols) {

	const int C = *std::max_element(p.color.begin(), p.color.end())+1;
	std::vector<bool> incomp(p.n());
	for(int u0 = 0; u0<int(p.n()); ++u0) {
		// int nedges = 0;
		std::vector<int> comp;
		std::set<int> neigh;
		std::vector<int> col_node(C, -1);
		comp.push_back(u0);
		incomp[u0] = 1;
		col_node[p.color[u0]] = u0;
		for(int v : g[u0]) {
			neigh.insert(v);
			// nedges += 1;
		}
		for(;;) {
			int minn = p.n()*2;
			int minu = -1;
			for(const auto u : neigh)
				if(col_node[p.color[u]]==-1) {
					int unn = 0;
					for(const auto v : g[u]) {
						if(incomp[v]) unn -= 1;
						else unn += 1;
					}
					if(unn<minn) {
						minn = unn;
						minu = u;
					}
				}
			if(minu==-1) break;
			// nedges += minn;
			col_node[p.color[minu]] = minu;
			comp.push_back(minu);
			incomp[minu] = 1;
			neigh.erase(minu);
			for(const auto v : g[minu])
				if(!incomp[v]) neigh.insert(v);

			std::sort(comp.begin(),comp.end());
			if(comp.size()>1) // maybe filter based on |inside edges|/|border edges| or something
				cols.push_back(comp);
		}
		for(const auto u : comp) incomp[u] = 0;
	}
}

void gen_with_karger_heuristic(const problem& p, std::vector<std::vector<int>>& cols) {
	auto c = heuristic_karger(p);
	for(auto &col: c) std::sort(col.begin(),col.end());
	std::move(c.begin(), c.end(), std::back_inserter(cols));
}

std::vector<std::vector<int>> gen_hotstart_columns(const problem &p) {
	std::vector<std::vector<int>> cols;
	std::vector<std::vector<int>> g(p.n());
	for(auto &[u,v]: p.edges) {
		g[u].push_back(v);
		g[v].push_back(u);
	}
	gen_with_greedy(p, g, cols);
	gen_with_karger_heuristic(p, cols);
	return cols;
}
