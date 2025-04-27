#pragma once

#include <cstdint>
#include <cassert>
#include <limits>
#include <vector>
#include <array>
#include <cstdint>

struct problem {
	std::vector<uint32_t> color;
	std::vector<std::array<uint32_t,2> > edges;
	size_t n() {
		return color.size();
	}
};

struct mcolumn {
    int value=0;
    std::vector<int> nodes;
};

struct Rng {
	Rng(uint64_t seed): state(seed) {}
	uint64_t state;
	uint64_t operator()() {
		state += 0x9e3779b97f4a7c15;
        uint64_t z = state;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
        z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
        return z ^ (z >> 31);
	}
	double f01() {
		return double((long double)(operator()())/(long double)(std::numeric_limits<uint64_t>::max()));
	}
};

using WT=uint64_t;
constexpr size_t WS = sizeof(WT)*8;
struct bitvec {
	size_t n;
	std::vector<WT> data;
	bitvec(const size_t N): n(N), data((N+WS-1)/WS,0) {}
	void zero() {
		std::fill(data.begin(),data.end(),0);
	}
	size_t count_ones() const {
		size_t ans=0;
		for(size_t i=0; i<data.size()-1; ++i) ans+=__builtin_popcountll(data[i]);
		ans+=__builtin_popcountll(data.back()&((1<<(n%WS))-1));
		return ans;
	}
	bitvec operator&(const bitvec& o) const {
		assert(n==o.n);
		bitvec ans(n);
		for(size_t i=0; i<data.size(); ++i) ans.data[i]=data[i]&o.data[i];
		return ans;
	}
	bitvec operator|(const bitvec& o) const {
		assert(n==o.n);
		bitvec ans(n);
		for(size_t i=0; i<data.size(); ++i) ans.data[i]=data[i]|o.data[i];
		return ans;
	}
	bitvec operator^(const bitvec& o) const {
		assert(n==o.n);
		bitvec ans(n);
		for(size_t i=0; i<data.size(); ++i) ans.data[i]=data[i]^o.data[i];
		return ans;
	}
	bitvec operator~() const {
		bitvec ans(n);
		for(size_t i=0; i<data.size(); ++i) ans.data[i]=~data[i];
		return ans;
	}
	void operator&=(const bitvec& o) {
		assert(n==o.n);
		for(size_t i=0; i<data.size(); ++i) data[i]&=o.data[i];
	}
	void operator|=(const bitvec& o) {
		assert(n==o.n);
		for(size_t i=0; i<data.size(); ++i) data[i]|=o.data[i];
	}
	void operator^=(const bitvec& o) {
		assert(n==o.n);
		for(size_t i=0; i<data.size(); ++i) data[i]^=o.data[i];
	}
	inline bool operator[](const size_t i) const {
		assert(i<n);
		return (data[i/WS]>>(i%WS))&1;
	}
	void set1(const size_t i) {
		assert(i<n);
		data[i/WS]|=(WT(1)<<(i%WS));
	}
	void set0(const size_t i) {
		assert(i<n);
		data[i/WS]&=~(WT(1)<<(i%WS));
	}
	void flip(const size_t i) {
		assert(i<n);
		data[i/WS]^=(WT(1)<<(i%WS));
	}
	void flipall() {
		for(size_t i=0; i<data.size(); ++i) data[i]=~data[i];
	}
	void set(const size_t i, bool b) {
		assert(i<n);
		data[i/WS]&=~(WT(1)<<(i%WS));
		data[i/WS]^=(WT(b)<<(i%WS));
	}
	void swap(const size_t i, const size_t j) {
		assert(i<n && j<n);
		bool bj = operator[](j);
		set(j,operator[](i));
		set(i,bj);
	}
	void pop_back() {
		--n;
		data.resize((n+WS-1)/WS);
	}
};

/// returns true if there's no element in common between a and b
bool compatible(const std::vector<int>& a, const std::vector<int>& b) {
    if(a.size()>b.size()) return compatible(b,a);
    auto bb = b.begin();
    for(auto &x: a) {
        bb = std::lower_bound(bb,b.end(),x);
        if(bb==b.end()) return 1;
        if(*bb==x) return 0;
    }
    return 1;
}


