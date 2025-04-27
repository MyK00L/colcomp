#pragma once

#include <lp_data/HConst.h>
#include <lp_data/HighsOptions.h>
#include <lp_data/HighsStatus.h>
#include <vector>
#include <Highs.h>
#include "mgraph.cpp"

class Master {
    public:
    mgraph g;
    std::vector<mcolumn> columns; // nodes in each column
    Highs highs; // lp
    vector<double> stab;
    double lagrange;
    Master(mgraph, std::vector<mcolumn>);
    // Master(const Master& o);
    void update_lp();
    void setup_highs();
    /// generates an initial set of columns
    void hot_start();
    /// tries to generate new columns
    /// returns true if new columns were generated
    bool generate_columns();
    /// keeps generating columns until optimality
    /// (TODO: or some condition is met, like lagrange passing the primal bound or being close to cg)
    double run();
};
