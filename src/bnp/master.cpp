#pragma once

#include <limits>
#include <lp_data/HConst.h>
#include <lp_data/HighsOptions.h>
#include <lp_data/HighsStatus.h>
#include <vector>
#include <iostream>
#include <Highs.h>
#include "bnp/master.hpp"

#include "bnp/pricing/heuristic_pricing_beamsearch.cpp"
#include "pricing/heuristic_pricing_annealing.cpp"
#include "pricing/exact_pricing_bnb.cpp"
#include "pricing/exact_pricing.cpp"

#include "pricing/pricing.hpp"

/*
idea: columns always have all colors, but master can choose which to actually use?
pricing idea: remove incrementally all nodes that have a multiplier bigger than the number of edges
pricing dual: then take the biggest degree-multiplier for each color
*/

/*
maximum number of partitions with all colors with the same number of nodes c:
c groups of c nodes -> >= 2c
min(n,2*maxc)??
*/


// do everything with SAT (?)

Master::Master(mgraph _mg, std::vector<mcolumn> _columns): g(_mg), columns(_columns), stab(g.n(),0), lagrange(std::numeric_limits<double>::infinity()) {}

void Master::setup_highs() {
    highs.clear();
    highs.setOptionValue("output_flag", false);
    // highs.setOptionValue("solver", "ipm");
    // highs.setOptionValue("seed",42); // TODO seed
    // highs.setMatrixFormat(MatrixFormat::kColwise);
    highs.changeObjectiveSense(ObjSense::kMaximize);

    vector<vector<int> > rows(g.n());
    for(int i=0; i<(int) columns.size(); ++i) {
        highs.addVar(0.0,1.0);
        highs.changeColCost(int(i), columns[i].value);
        for(auto &j:columns[i].nodes)rows[j].push_back(i);
    }
    vector<double> ones(columns.size(),1);
    for(size_t i=0; i<rows.size(); ++i) {
        highs.addRow(0,1,(int)rows[i].size(),rows[i].data(),ones.data());
    }
}

void Master::update_lp() {
    HighsStatus hs = highs.run();
    if(hs != HighsStatus::kOk) {
        std::cerr << "Highs status: "<<highsStatusToString(hs)<<std::endl;
    }
}

void Master::hot_start() {
    vector<double> ones(g.n(),1);
    auto cols = hot_start_cols(g);
    for(auto & col: cols) {
        columns.push_back({col.value, col.nodes});
        highs.addCol(col.value, 0, 1, (int)col.nodes.size(), col.nodes.data(), ones.data());
    }
}
// TODO: clean everything, better exact pricing
/// tries to generate new columns
/// returns true if new columns were generated
bool Master::generate_columns() {
    vector<double> duals = highs.getSolution().row_dual;
    
    bool is_exact=0;
    std::vector<mcolumn> pricing_columns;
    if(true) {
        pricing_problem pp(g,duals);
        pricing_columns = heuristic_pricing_beamsearch(pp);
        if(pricing_columns.empty()) pricing_columns = heuristic_pricing_annealing(pp);
    }
    if(pricing_columns.empty()) {
        for(int i=0; i<g.n(); ++i) { // single parameter stabilization
            duals[i]=duals[i]*0.8+stab[i]*0.2;
        }
        pricing_problem pp(g,duals);
        pricing_columns = exact_pricing_bnb(pp);
        // pricing_columns = exact_pricing(pp);
        /*
        auto check = exact_pricing(pp);
        std::cerr<<check[0].value<<' '<<pricing_columns[0].value<<'\n';
        */
        is_exact = 1;
    }

    double dual_sum=0;
    for(auto &i: duals) dual_sum+=i;

    vector<double> ones(g.n(),1);
    if(pricing_columns.size()) {
        for(auto& col: pricing_columns) {
            if(is_exact) {
                // update lagrangean dual TODO: check
                double tmp=col.value;
                for(auto &i: col.nodes) tmp-=duals[i];
                tmp*=g.n(); // should be * number of partitions (?)
                // tmp*=(g.n()+col.nodes.size()-1)/col.nodes.size(); // <- prob wrong, check
                // tmp*=double(g.n())/col.nodes.size();
                tmp+=dual_sum;
                // tmp+=g.value;
                // std::cerr<<"curr lagrange: "<<tmp<<'\n';
                if(tmp<lagrange) {
                    lagrange=tmp;
                    stab=duals;
                }
            }
            columns.push_back(col);
            highs.addCol((double)col.value, 0, 1, (int)col.nodes.size(), col.nodes.data(), ones.data());
        }
        return true;
    }

    return false;
}
/// keeps generating columns until optimality
/// (TODO: or some condition is met, like lagrange passing the primal bound or being close to cg)
double Master::run() {
    int it=0;
	do {
        update_lp();
        if(it%10==0) std::cerr<<it<<": "<<highs.getObjectiveValue()<<" | "<<lagrange<<std::endl;
        // if(int(lagrange) == int(g.value+highs.getObjectiveValue()+1e-6)) break;
        ++it;
    } while(generate_columns());
    std::cerr<<it<<": "<<highs.getObjectiveValue()<<" | "<<lagrange<<std::endl;
    // return highs.getObjectiveValue()+g.value;
	return lagrange;
}

// TODO: heuristic primal?
std::pair<int, std::vector<std::vector<int> > > Master::primal_ilp() {
    // set integrality constraints for columns
	for(int i=0; i<(int) columns.size(); ++i) highs.changeColIntegrality(i,HighsVarType::kInteger);
	highs.run();
	int obj = int(highs.getObjectiveValue()+0.1);
	const std::vector<double>& cv = highs.getSolution().col_value;
	std::vector<std::vector<int> > ans;
	for(int i=0; i<(int) columns.size(); ++i) if(cv[i]>0.5) {
		ans.push_back({});
		for(auto &u: columns[i].nodes) for(auto &ou: g.groups[u]) ans.back().push_back(ou);
	}

	// reset integrality constraints just in case
    for(int i=0; i<(int) columns.size(); ++i) highs.changeColIntegrality(i,HighsVarType::kContinuous);
	
	return {obj,ans};
}

std::vector<std::vector<int> > Master::get_columns() const {
	std::vector<std::vector<int> > ans;
	for(auto& [_v,col]: columns) {
		std::vector<int> ocol;
		for(auto &u:col) for(auto &ou: g.groups[u]) ocol.push_back(ou);
		std::sort(ocol.begin(), ocol.end());
		ans.push_back(ocol);
	}
	return ans;
}

