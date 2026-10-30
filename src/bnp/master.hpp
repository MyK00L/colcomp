#pragma once

#include "mgraph.cpp"
#include <Highs.h>
#include <lp_data/HConst.h>
#include <lp_data/HighsOptions.h>
#include <lp_data/HighsStatus.h>
#include <vector>

class Master {
  public:
	mgraph g;
	std::vector<mcolumn> columns; // nodes in each column
	Highs highs;				  // lp
	vector<double> stab;		  // single-parameter satbilization stability center
	double lagrange;			  // lagrangean dual value
	Master(mgraph, std::vector<mcolumn>);
	// Master(const Master& o);
	void update_lp();
	void setup_highs();
	/// tries to generate new columns
	/// returns true if new columns were generated
	bool generate_columns();
	/// keeps generating columns until a valid integer dual bound is found
	double run(int global_primal);
	/// returns the best partitioning with the current columns
	std::pair<int, std::vector<std::vector<int>>> primal_ilp(int primal);
	std::vector<std::vector<int>> get_columns() const;
};
