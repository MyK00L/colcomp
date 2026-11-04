#pragma once

#include "bnp/master.hpp"
#include <Highs.h>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <lp_data/HConst.h>
#include <vector>

#include "bnp/pricing/heuristic_pricing_beamsearch.cpp"
#include "pricing/exact_pricing_bnb.cpp"
#include "pricing/exact_pricing_bnb2.cpp"
#include "pricing/exact_pricing_ilp.cpp"

#include "pricing/pricing.hpp"

/*
MYSTERY:
If I set columns to be in [0,1] (where 1 is redundant),
Highs thinks it solved the rlmp optimally when it didn't (has positive col_duals),
and I generate duplicate columns :(
*/

enum ExactPricingMode {
	ILP,
	BNB,
	BNB_LPS,
} exact_pricing_mode = ExactPricingMode::BNB;

Master::Master(mgraph _mg, std::vector<mcolumn> _columns)
	: g(_mg), columns(_columns), stab(g.n(), 0), lagrange(std::numeric_limits<double>::infinity()) {}

/// setup the RLMP
void Master::setup_highs() {
	highs.clear();
	highs.setOptionValue("output_flag", false);
	highs.setOptionValue("solver", "simplex");
	highs.setOptionValue("simplex_strategy", "primal");
	highs.changeObjectiveSense(ObjSense::kMaximize);

	vector<vector<int>> rows(g.n());
	for(int i = 0; i<(int)columns.size(); ++i) {
		highs.addVar(0.0, kHighsInf);
		highs.changeColCost(int(i), columns[i].value);
		for(auto &j : columns[i].nodes) rows[j].push_back(i);
	}
	vector<double> ones(columns.size(), 1);
	for(size_t i = 0; i<rows.size(); ++i) {
		highs.addRow(-kHighsInf, 1, (int)rows[i].size(), rows[i].data(), ones.data());
	}
}
/// solve the RLMP
void Master::update_lp() {
	HighsStatus hs = highs.run();
	assert(highs.getModelStatus()==HighsModelStatus::kOptimal || highs.getModelStatus()==HighsModelStatus::kModelEmpty);
	assert(hs==HighsStatus::kOk);
	const auto &dc = highs.getSolution().col_dual;
	if(!dc.empty()) {
		double dcm = *std::max_element(dc.begin(), dc.end());
		if(dcm>=1e-5) std::cerr<<"WARNING: dcm="<<dcm<<std::endl;
		assert(dcm<1e-5);
	}
}

int pricing_steps = 0;
bool use_heuristic_pricing = 0;
/// tries to generate new columns
/// returns true if new columns were generated
bool Master::generate_columns() {
	pricing_steps += 1;
	vector<double> duals = highs.getSolution().row_dual;

	bool is_exact = 0;
	std::vector<mcolumn> pricing_columns;
	if(use_heuristic_pricing) { // heuristic pricing
		pricing_problem pp(g, duals);
		pricing_columns = heuristic_pricing_beamsearch(pp);
	}
	if(pricing_columns.empty()) {					   // exact pricing
		for(int i = 0; i<g.n(); ++i) {			   // single parameter stabilization
			// duals[i] = duals[i]*0.8+stab[i]*0.2; // value of 0.2 was found empirically
		}
		pricing_problem pp(g, duals);
		if(exact_pricing_mode==ExactPricingMode::ILP) {
			pricing_columns = exact_pricing_ilp(pp);
		} else if(exact_pricing_mode==ExactPricingMode::BNB) {
			pricing_columns = exact_pricing_bnb(pp);
		} else if(exact_pricing_mode==ExactPricingMode::BNB_LPS) {
			pricing_columns = exact_pricing_bnb2(pp);
		}
		is_exact = 1;
	}

	double dual_sum = 0; // sum(lambda)
	for(auto &i : duals) dual_sum += i;

	const vector<double> ones(g.n(), 1);
	if(pricing_columns.size()) {
		if(is_exact) {						// try updating lagrangean dual
			const int n_partitions = g.n(); // could use a better bound, but it would not matter
			// LR obj: max sum_{i}{lambda_i}+|P|*( sum_{ei}{e(ei)}-sum_{i}{lamda_i*x_i} )
			double tmp = pricing_columns[0].value; // first column has to be the optimal one
			for(auto &i : pricing_columns[0].nodes) tmp -= duals[i];
			tmp *= n_partitions;
			tmp += dual_sum;
			if(tmp<lagrange) { // found a new lagrangean dual bound
				lagrange = tmp;
				stab = duals; // update single parameter stabilization vector
			}
		}
		for(const auto &col : pricing_columns) { // add all generated columns to the rlmp
			columns.push_back(col); // TODO: check if they have positive rc when using single parameter stabilization
			highs.addCol((double)col.value, 0, kHighsInf, (int)col.nodes.size(), col.nodes.data(), ones.data());
		}

		return true;
	}

	return false;
}
/// keeps generating columns until a valid integer dual bound is found
double Master::run(int global_primal) {
	int it = 0;
	do {
		update_lp();
		assert(highs.getObjectiveValue()<=lagrange+1e-7);
		if(it%16==0) std::cerr<<it<<": "<<highs.getObjectiveValue()<<" ~ "<<lagrange<<std::endl;
		if(std::isfinite(lagrange)) {
			if(int(lagrange)==int(highs.getObjectiveValue()+1e-6)) {
				std::cerr<<it<<": "<<highs.getObjectiveValue()<<" ~ "<<lagrange<<" (exit1)"<<std::endl;
				return lagrange;
			}
			if(int(lagrange)<=global_primal) {
				std::cerr<<it<<": "<<highs.getObjectiveValue()<<" ~ "<<lagrange<<" (exit2)"<<std::endl;
				return lagrange;
			}
		}
		++it;
	} while(generate_columns());
	std::cerr<<it<<": "<<highs.getObjectiveValue()<<" ~ "<<lagrange<<" (exit3)"<<std::endl;
	// return lagrange;
	return highs.getObjectiveValue();
}
/// computes primal by solving integer restricted master problem for this master instance (not global primal, need to
/// add base)
std::pair<int, std::vector<std::vector<int>>> Master::primal_ilp(int primal) {
	// set integrality constraints for columns
	for(int i = 0; i<(int)columns.size(); ++i) highs.changeColIntegrality(i, HighsVarType::kInteger);
	{ // set objective>=primal+1
		std::vector<int> inds(columns.size());
		std::vector<double> vals(columns.size());
		for(int i = 0; i<(int)columns.size(); ++i) {
			inds[i] = i;
			vals[i] = columns[i].value;
		}
		highs.addRow(primal+1, kHighsInf, columns.size(), inds.data(), vals.data());
	}

	highs.run();
	assert(highs.getModelStatus()==HighsModelStatus::kOptimal ||
		   highs.getModelStatus()==HighsModelStatus::kInfeasible);
	int obj = highs.getModelStatus()==HighsModelStatus::kOptimal ? int(highs.getObjectiveValue()+0.1) : -1;
	const std::vector<double> &cv = highs.getSolution().col_value;
	std::vector<std::vector<int>> ans;
	vector<bool> taken(g.n());
	for(int i = 0; i<(int)columns.size(); ++i)
		if(cv[i]>0.5) {
			ans.push_back({});
			for(auto &u : columns[i].nodes) {
				taken[u] = 1;
				for(auto &ou : g.groups[u]) { ans.back().push_back(ou); }
			}
		}
	for(int u = 0; u<g.n(); ++u)
		if(!taken[u]) {
			ans.push_back({});
			for(auto &ou : g.groups[u]) { ans.back().push_back(ou); }
		}

	// reset integrality constraints and remove obj>=primal+1 constraint
	for(int i = 0; i<(int)columns.size(); ++i) highs.changeColIntegrality(i, HighsVarType::kContinuous);
	highs.deleteRows(highs.getNumRows()-1, highs.getNumRows()-1);

	return {obj, ans};
}
/// returns the columns generated by this master instance
std::vector<std::vector<int>> Master::get_columns() const {
	// TODO: throw away "bad" columns
	std::vector<std::vector<int>> ans;
	for(auto &[_v, col] : columns) {
		std::vector<int> ocol;
		for(auto &u : col)
			for(auto &ou : g.groups[u]) ocol.push_back(ou);
		std::sort(ocol.begin(), ocol.end());
		ans.push_back(ocol);
	}
	return ans;
}
