#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

struct problem {
	std::vector<uint32_t> color;
	std::vector<std::array<uint32_t, 2>> edges;
	size_t n() const { return color.size(); }

	void debug_print() {
		for(int u = 0; u<(int)n(); ++u) { std::cerr<<u<<'('<<color[u]<<")\n"; }
		for(const auto [u, v] : edges) { std::cerr<<u<<'('<<color[u]<<") "<<v<<'('<<color[v]<<")\n"; }
	}
};

struct mcolumn {
	int value = 0;
	std::vector<int> nodes;
};

struct Rng {
	Rng(uint64_t seed) : state(seed) {}
	uint64_t state;
	uint64_t operator()() {
		state += 0x9e3779b97f4a7c15;
		uint64_t z = state;
		z = (z^(z>>30))*0xbf58476d1ce4e5b9;
		z = (z^(z>>27))*0x94d049bb133111eb;
		return z^(z>>31);
	}
	double f01() { return double((long double)(operator()())/(long double)(std::numeric_limits<uint64_t>::max())); }
};

using WT = uint64_t;
constexpr size_t WS = sizeof(WT)*8;
struct bitvec {
	size_t n;
	std::vector<WT> data;
	bitvec(const size_t N) : n(N), data((N+WS-1)/WS, 0) {}
	void zero() { std::fill(data.begin(), data.end(), 0); }
	size_t count_ones() const {
		size_t ans = 0;
		for(size_t i = 0; i<data.size()-1; ++i) ans += __builtin_popcountll(data[i]);
		ans += __builtin_popcountll(data.back()&((uint64_t(1)<<(n%WS))-1));
		return ans;
	}
	bitvec operator&(const bitvec &o) const {
		assert(n==o.n);
		bitvec ans(n);
		for(size_t i = 0; i<data.size(); ++i) ans.data[i] = data[i]&o.data[i];
		return ans;
	}
	bitvec operator|(const bitvec &o) const {
		assert(n==o.n);
		bitvec ans(n);
		for(size_t i = 0; i<data.size(); ++i) ans.data[i] = data[i]|o.data[i];
		return ans;
	}
	bitvec operator^(const bitvec &o) const {
		assert(n==o.n);
		bitvec ans(n);
		for(size_t i = 0; i<data.size(); ++i) ans.data[i] = data[i]^o.data[i];
		return ans;
	}
	bitvec operator~() const {
		bitvec ans(n);
		for(size_t i = 0; i<data.size(); ++i) ans.data[i] = ~data[i];
		return ans;
	}
	void operator&=(const bitvec &o) {
		assert(n==o.n);
		for(size_t i = 0; i<data.size(); ++i) data[i] &= o.data[i];
	}
	void operator|=(const bitvec &o) {
		assert(n==o.n);
		for(size_t i = 0; i<data.size(); ++i) data[i] |= o.data[i];
	}
	void operator^=(const bitvec &o) {
		assert(n==o.n);
		for(size_t i = 0; i<data.size(); ++i) data[i] ^= o.data[i];
	}
	bool has_intersection(const bitvec &o) const {
		for(size_t i = 0; i<data.size()-1; ++i)
			if(data[i]&o.data[i]) return 1;
		if((data.back()&((uint64_t(1)<<(n%WS))-1))&(o.data.back()&((uint64_t(1)<<(n%WS))-1)))
			return 1;
		return 0;
	}
	inline bool operator[](const size_t i) const {
		assert(i<n);
		return (data[i/WS]>>(i%WS))&1;
	}
	void set1(const size_t i) {
		assert(i<n);
		data[i/WS] |= (WT(1)<<(i%WS));
	}
	void set0(const size_t i) {
		assert(i<n);
		data[i/WS] &= ~(WT(1)<<(i%WS));
	}
	void flip(const size_t i) {
		assert(i<n);
		data[i/WS] ^= (WT(1)<<(i%WS));
	}
	void flipall() {
		for(size_t i = 0; i<data.size(); ++i) data[i] = ~data[i];
	}
	void set(const size_t i, bool b) {
		assert(i<n);
		data[i/WS] &= ~(WT(1)<<(i%WS));
		data[i/WS] ^= (WT(b)<<(i%WS));
	}
	void swap(const size_t i, const size_t j) {
		assert(i<n && j<n);
		bool bj = operator[](j);
		set(j, operator[](i));
		set(i, bj);
	}
	void pop_back() {
		--n;
		data.resize((n+WS-1)/WS);
	}
	struct const_iterator {
		using iterator_category = std::forward_iterator_tag;
		const std::vector<WT> &data;
		const size_t n;
		uint64_t b;
		size_t idx;
		const_iterator(const std::vector<uint64_t> &_data, const size_t _n) : data(_data), n(_n), idx(0) {
			b = data[0];
			while(b==0 && idx<data.size()) {
				idx += 1;
				if(idx!=data.size()) b = data[idx];
			}
			if(idx!=data.size() && this->operator*()>=n) {
				idx = data.size();
				b = 0;
			}
		}
		const_iterator(const std::vector<uint64_t> &_data, const size_t _n, uint64_t _b, size_t _idx)
			: data(_data), n(_n), b(_b), idx(_idx) {}
		const_iterator operator++() {
			b ^= (b&(-b));
			while(b==0 && idx<data.size()) {
				idx += 1;
				if(idx!=data.size()) b = data[idx];
			}
			if(idx!=data.size() && this->operator*()>=n) {
				idx = data.size();
				b = 0;
			}
			return *this;
		}
		size_t operator*() { return idx*64+__builtin_ctzll(b&(-b)); }
		bool operator!=(const const_iterator &o) const { return idx!=o.idx || b!=o.b; }
	};
	const_iterator begin() const { return const_iterator(data, n); }
	const_iterator end() const { return const_iterator{data, n, 0, data.size()}; }
};

/// returns true if there's no element in common between a and b
bool compatible(const std::vector<int> &a, const std::vector<int> &b) {
	if(a.size()==1 && b.size()==1) return a[0]!=b[0]; // fastpath
	if(a.size()>b.size()) return compatible(b, a);
	auto bb = b.begin();
	for(auto &x : a) {
		bb = std::lower_bound(bb, b.end(), x);
		if(bb==b.end()) return 1;
		if(*bb==x) return 0;
	}
	return 1;
}
