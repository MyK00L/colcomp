#pragma once

#include "common.cpp"
#include "mgraph.cpp"
#include "pricing.hpp"
#include <Highs.h>
#include <vector>

/*
n_0 n_1 n_2 n_3...n_i
forall i: sum_j x[i,j]=n_i
forall j1,j2 max_i min(x[i,j1],x[i,j2])=1
max z = sum_j /\_i(x[i,j])

also doable with SAT, I guess
*/
std::vector<mcolumn> exact_pricing_ilp(const pricing_problem &p) {
	const double inf = kHighsInf;
	Highs h;
	h.setOptionValue("output_flag", false);
	h.changeObjectiveSense(ObjSense::kMaximize);
	int coli = 0;
	for(size_t i = 0; (int)i<p.g.n(); ++i) {
		h.addVar(0, 1); // x_i
		h.changeColIntegrality(coli, HighsVarType::kInteger);
		++coli;
	}
	for(int k = 0; k<p.g.C(); ++k) {
		vector<int> ind;
		vector<double> val;
		for(int i = 0; i<p.g.n(); ++i) {
			bool contains = 0;
			for(auto &ki : p.g.c[i])
				if(ki==k) {
					contains = 1;
					break;
				}
			if(contains) {
				ind.push_back(i);
				val.push_back(1);
			}
		}
		// forall k: sum_{u: k in c(u)}{x_u}<=1
		if(ind.size()) h.addRow(-inf, 1, ind.size(), ind.data(), val.data());
	}

	// edge constraints
	std::map<std::array<int, 2>, int> eis;
	for(int i = 0; i<int(p.g.n()); ++i) {
		for(auto [j, w] : p.g.g[i])
			if(j>i && w>0) {
				h.addVar(0, 1); // e_(u,v)
				eis[{i, j}] = coli;
				eis[{j, i}] = coli;
				h.changeColCost(coli, w);
				coli += 1;
			}
		for(int ci = 0; ci<p.g.C(); ++ci) {
			std::vector<int> inds = {i};
			std::vector<double> vals = {-1.0};
			for(auto [j, w] : p.g.g[i])
				if(w>0) {
					bool hasc = 0;
					for(auto col : p.g.c[j]) hasc |= col==ci;
					if(hasc) {
						inds.push_back(eis[{i, j}]);
						vals.push_back(1.0);
					}
				}
			// forall u,k: sum_{v:k in c(v)}{e_(u,v)}<=x_u
			if(inds.size()>1) h.addRow(-inf, 0, inds.size(), inds.data(), vals.data());
		}
	}

	for(size_t i = 0; i<p.cost.size(); ++i) h.changeColCost((int)i, -p.cost[i]);
	h.run();
	if(h.getObjectiveValue()<=0) return {};
	std::vector<int> nz;
	std::vector<double> cv = h.getSolution().col_value;
	for(size_t i = 0; i<p.cost.size(); ++i) {
		if(cv[i]>0.8) nz.push_back(i);
	}
	return {mcolumn{p.g.value_for_col(nz), nz}};
}
