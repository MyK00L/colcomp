#include <array>
#include <fstream>
#include <iostream>

#include "bnb.cpp"
#include "bnp/bnp.cpp"
#include "bnp/master.cpp"
#include "bnp/pricing/exact_pricing_bnb.cpp"
#include "common.cpp"
#include "heuristic/heuristic_karger.cpp"
#include "ilp.cpp"
#include "preprocess.cpp"
#include "scorer.cpp"
#include "heuristic/local_search.cpp"

using namespace std;

/*
 *	MCC: minimize number of components min cn
 *	MEC: maximize transitive closure sum(cs*(cs-1))
 *	MOP: minimize number of removed edges sum(cut(c))
 */

/*
 *2 colors -> bipartite matching
 *
 */

problem input(string filename) {
	ifstream fin(filename);
	vector<uint32_t> col;
	vector<array<uint32_t, 2>> edg;
	set<array<uint32_t, 2>> edgset;
	while(!fin.eof()) {
		char c = 0;
		fin>>c;
		if(fin.eof()) break;
		if(c=='v') {
			uint32_t a, b;
			fin>>a>>b;
			if(a+1>col.size()) col.resize(a+1);
			col[a] = b;
		} else if(c=='e') {
			uint32_t a, b;
			fin>>a>>b;
			if(a>b) swap(a, b);
			if(a>=col.size() || b>=col.size()) {
				std::cerr<<"node in edge has not been defined yet, ignoring"<<std::endl;
				continue;
			}
			if(col[a]==col[b]) {
				std::cerr<<"edge between same-color nodes in input, ignoring"<<std::endl;
				continue;
			}
			if(edgset.count({a, b})) {
				std::cerr<<"edge present twice in input, ignoring"<<std::endl;
				continue;
			}
			edgset.insert({a, b});
			edg.push_back({a, b});
		} else if(c=='#') {
			fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		} else {
			cerr<<"error reading input on character "<<int(c)<<" for file "<<filename<<endl;
			exit(1);
		}
	}
	return problem{col, edg};
}

int ccs(mgraph mg) {
	vector<bool> vis(mg.n(), 0);
	int res = 0;
	for(int i = 0; i<mg.n(); ++i)
		if(!vis[i]) {
			res += 1;
			queue<int> q;
			q.push(i);
			vis[i] = 1;
			while(!q.empty()) {
				int u = q.front();
				q.pop();
				for(auto &[v, _w] : mg.g[u])
					if(!vis[v]) {
						vis[v] = 1;
						q.push(v);
					}
			}
		}
	return res;
}
void reset_globals() {
	columns.clear();
	pricing_steps = 0;
	pricing_bnb_nodes = 0;
	pricing_bnb_runs = 0;
	pricing_bnb2_nodes = 0;
	pricing_bnb2_runs = 0;
	bnp_nodes = 0;
	use_heuristic_pricing = 0;
	exact_pricing_mode = ExactPricingMode::BNB;
}

void print_stats(const problem &p) {
	std::cout<<"N "<<p.n()<<std::endl;
	std::cout<<"M "<<p.edges.size()<<std::endl;
	std::cout<<"Qbound "<<max_useful_partitions(p)<<std::endl;
	mgraph mg(p);
	std::cout<<"C "<<mg.C()<<std::endl;
	auto ccs = mg.split_components();
	std::cout<<"ccs "<<ccs.size()<<std::endl;
	int maxn = 0, maxm = 0;
	for(auto &cc : ccs) {
		if(cc.n()>maxn) maxn = cc.n();
		if(cc.m()>maxm) maxm = cc.m();
	}
	std::cout<<"maxN "<<maxn<<std::endl;
	std::cout<<"maxM "<<maxm<<std::endl;
}

int main(int argc, char **argv) {
	cin.tie(0);
	ios_base::sync_with_stdio(0);
	if(argc<3) {
		std::string pname(argv[0]);
		std::cerr<<"Usage: "<<argv[0]
				 <<" [filename] [mode] (options)\n"
					 "possible modes:\n"
					 "stats : print graph info\n"
					 "merge : merge heuristic\n"
					 "karger : karger heuristic\n"
					 "ilp : solve MOP with ILP\n"
					 "bnp : solve MOP with BNP\n"
					 "preprocess : preprocess and stop\n"
					 "possible options:\n"
					 "nhs : no hotstart\n"
					 "nhp : no heuristic pricing\n"
					 "npp : no preprocessing before stats/merge/karger/ilp/bnp\n"
					 "ilp : use ilp for exact pricing\n"
					 "bnb : use default bnb for exact pricing\n"
					 "bnb_lps : use bnb with lp based dual bound for exact pricing"
				 <<std::endl;
		return 1;
	}
	string filename(argv[1]);
	problem p = input(filename);
	string mode(argv[2]);
	bool nhs = 0; // no hotstart
	bool npp = 0; // no preprocess
	use_heuristic_pricing = 1;
	for(int i = 3; i<argc; ++i) {
		std::string arg(argv[i]);
		if(arg=="nhs") nhs = 1;
		else if(arg=="nhp") use_heuristic_pricing = 0;
		else if(arg=="npp") npp = 1;
		else if(arg=="ilp") exact_pricing_mode = ExactPricingMode::ILP;
		else if(arg=="bnb") exact_pricing_mode = ExactPricingMode::BNB;
		else if(arg=="bnb_lps") exact_pricing_mode = ExactPricingMode::BNB_LPS;
		else {
			std::cerr<<"unknown option "<<arg<<std::endl;
			exit(1);
		}
	}
	clock_t t0 = clock();
	if(mode=="stats") {
		if(!npp) preprocess(p);
		print_stats(p);
	} else if(mode=="merge") {
		if(!npp) preprocess(p);
		auto sol = heuristic_merge(p);
		std::cout<<"score "<<score(p, sol)<<std::endl;
	} else if(mode=="karger") {
		if(!npp) preprocess(p);
		auto sol = heuristic_karger(p);
		std::cout<<"score "<<score(p, sol)<<std::endl;
	} else if(mode=="ilp") {
		if(!npp) preprocess(p);
		auto sol = ilp(p);
		std::cout<<"score "<<score(p, sol)<<std::endl;
	} else if(mode=="bnp") {
		if(!npp) preprocess(p);
		auto sol = bnp_main(p, !nhs);
		std::cout<<"score "<<score(p, sol)<<std::endl;
		std::cout<<"bnp_nodes "<<bnp_nodes<<std::endl;
		std::cout<<"pricing_steps "<<pricing_steps<<std::endl;
		std::cout<<"pricing_bnb_nodes "<<pricing_bnb_nodes<<std::endl;
		std::cout<<"pricing_bnb_runs "<<pricing_bnb_runs<<std::endl;
	} else if(mode=="preprocess") {
		preprocess(p);
	} else {
		std::cerr<<"unknown mode "<<mode<<std::endl;
		return 1;
	}
	clock_t t1 = clock();
	std::cout<<"time_millis "<<1000*uint64_t(t1-t0)/CLOCKS_PER_SEC<<std::endl;

	return 0;
}

// heuristic idea: keep cutting in the min-cut between two nodes with the same color
