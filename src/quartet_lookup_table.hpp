#pragma once

//#define USE_STXXL

#include <array>
#include <cassert>
#include <cstdint>
#include <vector>

#ifdef USE_STXXL
#include <stxxl/vector>
#endif

// =================================================================================================
//     Quartet Lookup Table
// =================================================================================================

template<typename LookupIntType>
class QuartetLookupTable {
public:

	// -------------------------------------------------------------------------
	//     Typedefs and Enums
	// -------------------------------------------------------------------------

	using QuartetTuple = std::array< LookupIntType, 3 >;

	// -------------------------------------------------------------------------
	//     Constructors and Rule of Five
	// -------------------------------------------------------------------------

	QuartetLookupTable() :
			num_taxa_(0) {
	}

	QuartetLookupTable(size_t num_taxa) {
		init(num_taxa);
	}

	~QuartetLookupTable() {
#ifdef USE_STXXL
		while (!quartet_lookup_.empty()) {
			quartet_lookup_.pop_back();
		}
#endif
	}

	QuartetLookupTable(QuartetLookupTable const&) = default;
	QuartetLookupTable(QuartetLookupTable&&) = default;

	QuartetLookupTable& operator=(QuartetLookupTable const&) = default;
	QuartetLookupTable& operator=(QuartetLookupTable&&) = default;

	// -------------------------------------------------------------------------
	//     Public Interface
	// -------------------------------------------------------------------------

	void init(size_t num_taxa) {
		num_taxa_ = num_taxa;
		// init_binom_lookup_(num_taxa);
		init_quartet_lookup_(num_taxa);
	}

	size_t num_taxa() const {
		return num_taxa_;
	}

	size_t size() const {
		return (quartet_lookup_.size() * 3 * sizeof(LookupIntType)) + (binom_lookup_.size() + 1) * sizeof(size_t);
	}

	QuartetTuple& get_tuple(size_t a, size_t b, size_t c, size_t d) {
		size_t const id = lookup_index_(a, b, c, d);
		assert(id < quartet_lookup_.size());
		return quartet_lookup_[id];
	}

	QuartetTuple const& get_tuple(size_t a, size_t b, size_t c, size_t d) const {
		size_t const id = lookup_index_(a, b, c, d);
		assert(id < quartet_lookup_.size());
		return quartet_lookup_[id];
	}

	size_t tuple_index(size_t a, size_t b, size_t c, size_t d) const {
		// Get all comparisons that we need.
		bool const ac = (a<c);
		bool const ad = (a<d);
		bool const bc = (b<c);
		bool const bd = (b<d);

		// Check first and third case. Second one is implied.
		bool const x = ((ac) & (ad) & (bc) & (bd)) | ((!ac) & (!bc) & (!ad) & (!bd));
		bool const ab_in_cd = ((!ac) & (ad) & (!bc) & (bd)) | ((!ad) & (ac) & (!bd) & (bc));
		bool const cd_in_ab = ((ac) & (!bc) & (ad) & (!bd)) | ((bc) & (!ac) & (bd) & (!ad));
		bool const z = ab_in_cd | cd_in_ab;
		bool const y = !x & !z;

		// Only one can be set.
		assert(!(x & y & z));
		assert(x ^ y ^ z);
		assert(x | y | z);
		size_t const r = static_cast<size_t>(y) + 2 * static_cast<size_t>(z);

		// Result has to be fitting.
		assert(r < 3);
		assert((x && !y && !z && r == 0) || (!x && y && !z && r == 1) || (!x && !y && z && r == 2));
		return r;
	}

	// -------------------------------------------------------------------------
	//     Private Members
	// -------------------------------------------------------------------------

	void init_binom_lookup_(size_t num_taxa) {
		binom_lookup_ = std::vector<size_t>(num_taxa * 5, 0);

		for (size_t i = 0; i < num_taxa; ++i) {
			for (size_t j = 0; j <= 4; ++j) {
				if (i == j || j == 0 || i == 0) {
					if (i == 0 && j > 0) {
						binom_lookup_[i * 5 + j] = 0;
					} else {
						binom_lookup_[i * 5 + j] = 1;
					}
				} else {
					binom_lookup_[i * 5 + j] = binom_lookup_[(i - 1) * 5 + j - 1] + binom_lookup_[(i - 1) * 5 + j];
				}
			}
		}
	}

	void init_quartet_lookup_(size_t num_taxa) {
		// calculate ncr(n, 4)
		size_t const n = (num_taxa * (num_taxa - 1) * (num_taxa - 2) * (num_taxa - 3)) / 24;
		quartet_lookup_ = std::vector<QuartetTuple>(n, {{ 0, 0, 0 }});
	}

	size_t binom_coefficient_sum_(size_t a, size_t b, size_t c, size_t d) const {
		// auto binom = [&]( size_t n, size_t k ) {
		// 	if( n * 5 + k >= binom_lookup_.size() ) {
		// 		std::cout << "n " << n << " k " << k << " size " << binom_lookup_.size() << "\n";
		// 	}
		//
		// 	assert( n * 5 + k < binom_lookup_.size() );
		// 	return binom_lookup_[ n * 5 + k ];
		// };
		//
		// // We expect sorted input, starting at the largest.
		// assert(a > b && b > c && c > d);
		//
		// assert(d == binom(d, 1));
		// size_t res = d;
		// res += binom(c, 2);
		// res += binom(b, 3);
		// res += binom(a, 4);

		// Alternative, without lookup.
		size_t res = 0;
		res += ( a * (a - 1) * (a - 2) * (a - 3) ) / 24;
		res += ( b * (b - 1) * (b - 2) ) / 6;
		res += ( c * (c - 1) ) / 2;
		res += d;

		return res;
	}

	size_t lookup_index_(size_t a, size_t b, size_t c, size_t d) const {
		size_t ta, tb, tc, td; // from largest to smallest
		size_t low1, high1, low2, high2, middle1, middle2;

		if (a < b) {
			low1 = a;
			high1 = b;
		} else {
			low1 = b;
			high1 = a;
		}
		if (c < d) {
			low2 = c;
			high2 = d;
		} else {
			low2 = d;
			high2 = c;
		}

		if (low1 < low2) {
			td = low1;
			middle1 = low2;
		} else {
			td = low2;
			middle1 = low1;
		}
		if (high1 > high2) {
			ta = high1;
			middle2 = high2;
		} else {
			ta = high2;
			middle2 = high1;
		}
		if (middle1 < middle2) {
			tc = middle1;
			tb = middle2;
		} else {
			tc = middle2;
			tb = middle1;
		}

		return binom_coefficient_sum_(ta, tb, tc, td);
	}

	// -------------------------------------------------------------------------
	//     Data Members
	// -------------------------------------------------------------------------

#ifdef USE_STXXL
	stxxl::VECTOR_GENERATOR<QuartetTuple>::result quartet_lookup_;
#else
	std::vector<QuartetTuple> quartet_lookup_;
#endif

	std::vector<size_t> binom_lookup_;

	size_t num_taxa_;

};
