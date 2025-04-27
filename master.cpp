#pragma once

#include <limits>
#include <lp_data/HConst.h>
#include <lp_data/HighsOptions.h>
#include <lp_data/HighsStatus.h>
#include <vector>
#include <iostream>
#include <Highs.h>
#include "master.hpp"

#include "heuristic_pricing_beamsearch.cpp"
#include "heuristic_pricing_annealing.cpp"
#include "exact_pricing_bnb.cpp"

#include "pricing.hpp"

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

Master::Master(mgraph _mg, std::vector<mcolumn> _columns): g(_mg), columns(_columns), lagrange(std::numeric_limits<double>::infinity()) {
    // setup_highs();
}
//Master::Master(const Master&o):g(o.g),columns(o.columns),lagrange(std::numeric_limits<double>::infinity()) {
//    setup_highs();
//}
void Master::setup_highs() {
    highs.clear();
    highs.setOptionValue("solver", "ipm");
    highs.setOptionValue("output_flag", false);
    //highs.setOptionValue("seed",42); // TODO seed
    highs.setMatrixFormat(MatrixFormat::kColwise);
    highs.changeObjectiveSense(ObjSense::kMaximize);
    vector<vector<int> > rows(g.n());
    for(int i=0; i<(int) columns.size(); ++i) {
        highs.addVar(0,1);
        highs.changeColCost(int(i), columns[i].value);
        for(auto &j:columns[i].nodes)rows[j].push_back(i);
    }
    vector<double> ones(columns.size(),1);
    for(size_t i=0; i<rows.size(); ++i) {
        highs.addRow(0,1,rows[i].size(),rows[i].data(),ones.data());
    }
}

void Master::update_lp() {
    HighsStatus hs = highs.run();
    if(hs != HighsStatus::kOk) {
        std::cerr << "Highs status: "<<highsStatusToString(hs)<<std::endl;
    }
}
/// tries to generate new columns
/// returns true if new columns were generated
// double nann = 0;
// double nbs = 0;
vector<double> stab;
bool phase2=0;
bool prephase2=0;

void Master::hot_start() {
    vector<double> duals(g.n(),0);
    vector<double> ones(g.n(),1);
    pricing_problem pp(g,duals);
    auto cols = hot_start_cols(pp);
    for(auto & col: cols) {
        columns.push_back({(int(col.value)), col.nz});
        highs.addCol(col.value, 0, 1, (int)col.nz.size(), col.nz.data(), ones.data());
    }
}
// TODO: clean everything, better exact pricing
bool Master::generate_columns() {
    vector<double> duals = highs.getSolution().row_dual;
    
    bool is_exact=0;
    std::vector<column> pricing_columns;
    if(false) {
        pricing_problem pp(g,duals);
        pricing_columns = heuristic_pricing_beamsearch(pp);
        if(pricing_columns.empty()) pricing_columns = heuristic_pricing_annealing(pp);
    }
    if(pricing_columns.empty()) {
        for(int i=0; i<g.n(); ++i) { // single parameter stabilization
            duals[i]=duals[i]*0.6+stab[i]*0.4;
        }
        pricing_problem pp(g,duals);
        pricing_columns = exact_pricing_bnb(pp);
        /*
        auto check = exact_pricing(pp);
        std::cerr<<check[0].value<<' '<<pricing_columns[0].value<<'\n';
        */
        is_exact = 1;
    } else prephase2=0;


    double dual_sum=0;
    for(auto &i: duals) dual_sum+=i;

    vector<double> ones(g.n(),1);
    if(pricing_columns.size()) {
        for(auto& col: pricing_columns) {
            double value = col.value;
            for(auto &i: col.nz) value+=duals[i];

            if(is_exact) {
                double tmp = dual_sum+col.value*g.n();//(dual_sum+col.value)*g.n(); // ?????
                // std::cerr<<"curr lagrange: "<<tmp<<'\n';
                if(tmp<lagrange) {
                    lagrange=tmp;
                    stab=duals;
                }
            }
            //std::cerr<<"cv2: "<<cv2<<'\n';
            //std::cerr<<"lagrange: "<<lagrange<<'\n';
            //std::cerr<<"dual_sum: "<<dual_sum<<'\n';
            //std::cerr<<"col val: "<<col.value<<"\ncol: ";
            //for(auto &i: col.nz) std::cerr<<i<<' ';
            //std::cerr<<'\n';
            //for(auto &i: col.nz) std::cerr<<duals[i]<<' ';
            //std::cerr<<std::endl;
            //std::cerr<<"aoeu "<<col.value<<' '<<value<<'\n';

            columns.push_back({int(round(value)), col.nz});
            highs.addCol(value, 0, 1, (int)col.nz.size(), col.nz.data(), ones.data());
        }
        return true;
    }
    
    return false;
}
/// keeps generating columns until optimality
/// (TODO: or some condition is met, like lagrange passing the primal bound or being close to cg)
double Master::run() {
    phase2 = prephase2 = 0;
    int it=0;
    stab=vector<double>(g.n(),0);
    do {
        update_lp();
        //std::cerr<<"objective: "<<g.value + highs.getInfo().objective_function_value<<std::endl;
        /*if(it%10==0) */std::cerr<<it<<": "<<g.value+highs.getObjectiveValue()<<" | "<<lagrange<<std::endl;
        if(int(lagrange) == int(g.value+highs.getObjectiveValue()+1e-6)) break;
        ++it;
    } while(generate_columns());
    // std::cerr<<"nbs:"<<nbs<<" nann:"<<nann<<std::endl;
    return highs.getObjectiveValue()+g.value;
}
