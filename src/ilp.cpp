#include "common.cpp"
#include "scorer.cpp"
#include <Highs.h>

int max_useful_partitions(const problem &p) {
	// compute bound on max useful number of partitions
	vector<int> nc;
	for(auto &c : p.color) {
		if(c>=nc.size()) nc.resize(c+1, 0);
		nc[c] += 1;
	}
	int max_edges = 0;
	for(auto &c : nc) { max_edges += c*(c-1)/2; }
	return std::min((int)p.n(), int(sqrt(max_edges*2+0.25)+0.5+0.001));
}

std::vector<std::vector<int>> ilp(problem p) {
	const double inf = kHighsInf;
	int n = (int)p.n();
	int P = max_useful_partitions(p);
	int m = (int)p.edges.size();
	int C = 1+*std::max_element(p.color.begin(), p.color.end());
	int nvar_part = n+m;
	Highs highs;
	highs.setOptionValue("output_flag", false);
	highs.changeObjectiveSense(ObjSense::kMaximize);
	// add variables
	for(int parti = 0; parti<P; ++parti) {
		for(int ni = 0; ni<n; ++ni) {
			highs.addVar(0, 1);
			highs.changeColIntegrality(parti*nvar_part+ni, HighsVarType::kInteger);
		}
		for(int ei = 0; ei<m; ++ei) {
			highs.addVar(0, 1);
			highs.changeColCost(nvar_part*parti+n+ei, 1);
		}
	}
	// add constraints
	for(int parti = 0; parti<P; ++parti) {
		/*
		for(int i=0; i<n; ++i) for(int j=i+1; j<n; ++j) if(p.color[i]==p.color[j]) { // color constraints
			int ind[2] = {parti*nvar_part+i,parti*nvar_part+j};
			double val[2] = {1,1};
			highs.addRow(-inf,1,2,ind,val); // (-inf <=) i+j<=1
		}
		*/
		for(int c = 0; c<C; ++c) { // color constraints
			vector<int> ind;
			vector<double> val;
			for(int i = 0; i<n; ++i)
				if((int)p.color[i]==c) {
					ind.push_back(parti*nvar_part+i);
					val.push_back(1);
				}
			if(ind.size()<=1) continue;
			highs.addRow(-inf, 1, ind.size(), ind.data(), val.data()); // (-inf <=) i+j+k+...<=1
		}
		for(int ei = 0; ei<m; ++ei) { // edge objective constraints
			int i = p.edges[ei][0];
			int j = p.edges[ei][1];
			int nz[2] = {parti*nvar_part+i, nvar_part*parti+n+ei};
			double val[2] = {-1, 1};
			highs.addRow(-inf, 0, 2, nz, val); // (-inf <=) e-i<=0
			nz[0] = parti*nvar_part+j;
			highs.addRow(-inf, 0, 2, nz, val); // (-inf <=) e-j<=0
		}
	}
	for(int i = 0; i<n; ++i) { // partition constraints
		vector<int> ind(P);
		vector<double> val(P, 1);
		for(int parti = 0; parti<P; ++parti) ind[parti] = parti*nvar_part+i;
		highs.addRow(-inf, 1, P, ind.data(), val.data()); // (-inf <=) i+i+...+i<=1
	}
	highs.run();
	std::cerr<<"ilp solution: "<<highs.getObjectiveValue()<<std::endl;

	std::vector<std::vector<int>> ans(P);
	{
		auto &cv = highs.getSolution().col_value;
		for(int parti = 0; parti<P; ++parti) {
			for(int ni = 0; ni<n; ++ni) {
				if(cv[parti*nvar_part+ni]>0.5) { ans[parti].push_back(ni); }
			}
		}
		for(int i = 0; i<(int)ans.size(); ++i) {
			if(ans[i].empty()) {
				ans[i].swap(ans.back());
				ans.pop_back();
				--i;
			}
		}
	}
	assert(int(highs.getObjectiveValue()+0.1)==score(p, ans));
	return ans;
}