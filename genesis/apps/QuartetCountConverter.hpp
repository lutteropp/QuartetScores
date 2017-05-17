#pragma once
#include <cmath>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <iostream>

/**
 * A quartet is a set {a,b,c,d} of four taxa.
 * This class stores the counts of observing the topologies ab|cd, ac|bd, and ad|bc in the evaluation tree set.
 * Note that the taxon IDs a, b, c, and d are not stored in this data structure.
 */
class QuartetCount {
public:
	QuartetCount() {
		count_12_34 = 0;
		count_13_24 = 0;
		count_14_23 = 0;
	}
	/**
	 * Get the count of the topology ab|cd in the evaluation tree set.
	 * @param a the first taxon ID
	 * @param b the second taxon ID
	 * @param c the third taxon ID
	 * @param d the fourth taxon ID
	 */
	size_t getCount(int a, int b, int c, int d) {
		if (a <= b && b <= c && c <= d) { // 12_34
			return count_12_34;
		} else if (a <= b && b <= d && d <= c) { // 12_43
			return count_12_34;
		} else if (b <= a && a <= c && c <= d) { // 21_34
			return count_12_34;
		} else if (b <= a && a <= d && d <= c) { // 21_43
			return count_12_34;
		} else if (c <= d && d <= a && a <= b) { // 34_12
			return count_12_34;
		} else if (c <= d && d <= b && b <= a) { // 34_21
			return count_12_34;
		} else if (d <= c && c <= a && a <= b) { // 43_12
			return count_12_34;
		} else if (d <= c && c <= b && b <= a) { // 43_21
			return count_12_34;

		} else if (a <= c && c <= b && b <= d) { // 13_24
			return count_13_24;
		} else if (a <= c && c <= d && d <= b) { // 13_42
			return count_13_24;
		} else if (c <= a && a <= b && b <= d) { // 31_24
			return count_13_24;
		} else if (c <= a && a <= d && d <= b) { // 31_42
			return count_13_24;
		} else if (b <= d && d <= a && a <= c) { // 24_13
			return count_13_24;
		} else if (d <= b && b <= a && a <= c) { // 42_13
			return count_13_24;
		} else if (b <= d && d <= c && c <= a) { // 24_31
			return count_13_24;
		} else if (d <= b && b <= c && c <= a) { // 42_31
			return count_13_24;

		} else {
			return count_14_23;
		}
	}

	/**
	 * Increase the count of the topology ab|cd by a given amount.
	 * @param a the first taxon ID
	 * @param b the second taxon ID
	 * @param c the third taxon ID
	 * @param d the fourth taxon ID
	 * @param count the amount to add to the count of ab|cd
	 */
	void addCount(int a, int b, int c, int d, size_t count) { // increase the count for ab|cd
		if (a <= b && b <= c && c <= d) { // 12_34
			count_12_34 += count;
		} else if (a <= b && b <= d && d <= c) { // 12_43
			count_12_34 += count;
		} else if (b <= a && a <= c && c <= d) { // 21_34
			count_12_34 += count;
		} else if (b <= a && a <= d && d <= c) { // 21_43
			count_12_34 += count;
		} else if (c <= d && d <= a && a <= b) { // 34_12
			count_12_34 += count;
		} else if (c <= d && d <= b && b <= a) { // 34_21
			count_12_34 += count;
		} else if (d <= c && c <= a && a <= b) { // 43_12
			count_12_34 += count;
		} else if (d <= c && c <= b && b <= a) { // 43_21
			count_12_34 += count;

		} else if (a <= c && c <= b && b <= d) { // 13_24
			count_13_24 += count;
		} else if (a <= c && c <= d && d <= b) { // 13_42
			count_13_24 += count;
		} else if (c <= a && a <= b && b <= d) { // 31_24
			count_13_24 += count;
		} else if (c <= a && a <= d && d <= b) { // 31_42
			count_13_24 += count;
		} else if (b <= d && d <= a && a <= c) { // 24_13
			count_13_24 += count;
		} else if (d <= b && b <= a && a <= c) { // 42_13
			count_13_24 += count;
		} else if (b <= d && d <= c && c <= a) { // 24_31
			count_13_24 += count;
		} else if (d <= b && b <= c && c <= a) { // 42_31
			count_13_24 += count;

		} else {
			count_14_23 += count;
		}
	}
private:
	size_t count_12_34; /**< count of ab|cd with a<b<c<d */
	size_t count_13_24; /**< count of ac|bd with a<b<c<d */
	size_t count_14_23; /**< count of ad|bc with a<c<b<d */
};

/**
 * This class maps a quartet {a,b,c,d} to a number and vice-versa.
 * We obtain a bijective mapping of all quartets to {0,...,N-1}, where N is the total number of quartets.
 */
class QuartetConverter {
public:
	size_t quartetToNumber(int aIdx, int bIdx, int cIdx, int dIdx);
	void numberToQuartet(size_t z, int &aIdx, int &bIdx, int &cIdx, int &dIdx);
	QuartetConverter(int numberOfTaxa);
	void test();
private:
	size_t f2(int a);
	size_t f3(int b);
	size_t f4(int c);
	size_t binarySearch(std::vector<size_t> &array, size_t z);

	int n; /**< number of taxa in the reference tree */
	std::vector<size_t> prefixSumF2; /**< prefixSumF2[i] := f2(0) + f2(1) + ... + f2(i) */
	std::vector<size_t> prefixSumF3; /**< prefixSumF3[i] := f3(0) + f3(1) + ... + f3(i) */
	std::vector<size_t> prefixSumF4; /**< prefixSumF4[i] := f4(0) + f4(1) + ... + f4(i) */

	std::vector<std::vector<size_t> > binom; /**< binomial coefficients */
};

/**
 * @param numberOfTaxa number of taxa in the reference tree
 */
QuartetConverter::QuartetConverter(int numberOfTaxa) {
	n = numberOfTaxa;
	// preprocess prefix sums
	prefixSumF2.resize(n - 2);
	prefixSumF3.resize(n - 2);
	prefixSumF4.resize(n - 2);
	prefixSumF2[0] = 1;
	prefixSumF3[0] = 1;
	prefixSumF4[0] = 1;
	for (int i = 1; i < n - 3; ++i) {
		prefixSumF2[i] = prefixSumF2[i - 1] + f2(i);
		prefixSumF3[i] = prefixSumF3[i - 1] + f3(i + 1);
		prefixSumF4[i] = prefixSumF4[i - 1] + f4(i + 2);
	}

	// preprocess binomial coefficients
	binom.resize(numberOfTaxa);
	for (int i = 0; i < numberOfTaxa; ++i) {
		binom[i].resize(5);
		for (int j = 0; j <= 4; ++j) {
			if (i == j || j == 0 || i == 0) {
				binom[i][j] = 0;
			} else {
				binom[i][j] = binom[i - 1][j - 1] + binom[i - 1][j];
			}
		}
	}

	// TODO: This is just for debug, remove this again.
	int a, b, c, d;
	numberToQuartet(1234, a, b, c, d);
	if (quartetToNumber(a, b, c, d) == 1234) {
		std::cout << "MAPPING TEST PASSED.\n";
	} else {
		std::cerr << "MAPPING TEST FAILED. Got " << quartetToNumber(a, b, c, d) << " instead of 1234.\n";
	}
}

/**
 * Convert the quartet {aIdx,bIdx,cIdx,dIdx} to a number, assuming lexicographical order of the DECREASING sequence.
 * @param aIdx a taxon ID
 * @param bIdx a taxon ID
 * @param cIdx a taxon ID
 * @param dIdx a taxon ID
 */
size_t QuartetConverter::quartetToNumber(int aIdx, int bIdx, int cIdx, int dIdx) {
	// TODO: This currently assumes lexicographical order of the DECREASING sequence.
	int a, b, c, d; // from largest to smallest
	int low1, high1, low2, high2, middle1, middle2;
	if (aIdx < bIdx) {
		low1 = aIdx;
		high1 = bIdx;
	} else {
		low1 = bIdx;
		high1 = aIdx;
	}
	if (cIdx < dIdx) {
		low2 = cIdx;
		high2 = dIdx;
	} else {
		low2 = dIdx;
		high2 = cIdx;
	}
	if (low1 < low2) {
		d = low1;
		middle1 = low2;
	} else {
		d = low2;
		middle1 = low1;
	}
	if (high1 > high2) {
		a = high1;
		middle2 = high2;
	} else {
		a = high2;
		middle2 = high1;
	}
	if (middle1 < middle2) {
		c = middle1;
		b = middle2;
	} else {
		c = middle2;
		b = middle1;
	}

	size_t res = 0;

	res += binom[a][4];
	res += binom[b][3];
	res += binom[c][2];
	res += binom[d][1];

	/*res += ( a * (a - 1) * (a - 2) * (a - 3) ) / 24;
	 res += ( b * (b - 1) * (b - 2) ) / 6;
	 res += ( c * (c - 1) ) / 2;
	 res += d;*/
	return res;
}

/**
 * Convert a number z to a quartet {aIdx, bIdx, cIdx, dIdx}, with the property that aIdx<bIdx<cIdx<dIdx.
 * This currently assumes lexicographical order of the INCREASING sequence.
 * @param z the number representing a quartet
 * @param aIdx a taxon ID
 * @param bIdx a taxon ID
 * @param cIdx a taxon ID
 * @param dIdx a taxon ID
 */
void QuartetConverter::numberToQuartet(size_t z, int &aIdx, int &bIdx, int &cIdx, int &dIdx) {
	// TODO: This currently assumes lexicographical order of the INCREASING sequence.
	z++; // because of start by 0 vs. start by 1

	size_t wantedT1 = z;
	aIdx = binarySearch(prefixSumF2, z) + 1;
	size_t foundT1 = prefixSumF2[aIdx - 1];

	if (wantedT1 == foundT1) {
		bIdx = aIdx + 1;
		cIdx = aIdx + 2;
		dIdx = aIdx + 3;
	} else {
		size_t wantedT2 = (prefixSumF3[aIdx - 1]) + (wantedT1 - foundT1);
		bIdx = binarySearch(prefixSumF3, wantedT2) + 2;
		size_t foundT2 = prefixSumF3[bIdx - 2];

		if (wantedT2 == foundT2) {
			cIdx = bIdx + 1;
			dIdx = bIdx + 2;
		} else {
			size_t wantedT3 = (prefixSumF4[bIdx - 2]) + (wantedT2 - foundT2);
			cIdx = binarySearch(prefixSumF4, wantedT3) + 3;
			size_t foundT3 = prefixSumF4[cIdx - 3];

			if (wantedT3 == foundT3) {
				dIdx = cIdx + 1;
			} else {
				dIdx = wantedT3 - foundT3 + cIdx + 1;
			}
		}
	}

	// because of start by 0 vs. start by 1
	aIdx--;
	bIdx--;
	cIdx--;
	dIdx--;
}

/**
 * Finds the highest position in the given array such that the entry is smaller than z.
 * @param array the array to search in
 * @param z the number to compare to
 */
size_t QuartetConverter::binarySearch(std::vector<size_t> &array, size_t z) {
	size_t first = 0;
	size_t last = n - 3;
	size_t middle = (first + last) / 2;
	size_t lastSmaller = 0;
	while (first < last) {
		if (array[middle] < z) {
			first = middle + 1;
			lastSmaller = middle;
		} else if (array[middle] > z) {
			last = middle - 1;
		} else { // array[middle] == z
			lastSmaller = middle;
			break;
		}
		middle = (first + last) / 2;
	}
	return lastSmaller;
}

size_t QuartetConverter::f2(int a) {
	long double nDouble = n;
	long double aDouble = a;
	long double res = (nDouble - aDouble) * (nDouble - 1 - aDouble) * (nDouble - 2 - aDouble) / 6;
	return round(res);
}

size_t QuartetConverter::f3(int b) {
	long double nDouble = n;
	long double bDouble = b;
	long double res = (nDouble - bDouble) * (nDouble - 1 - bDouble) / 2;
	return round(res);
}

size_t QuartetConverter::f4(int c) {
	return n - c;
}

void QuartetConverter::test() {
	// Quartet converter test
	int a, b, c, d;

	for (int i = 0; i <= 10; ++i) {
		numberToQuartet(i, a, b, c, d);
		std::cout << "Number " << i << " has quartet (" << a << "," << b << "|" << c << "," << d << ")" << std::endl;
	}

	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			for (int k = j + 1; k < n; ++k) {
				for (int l = k + 1; l < n; ++l) {
					size_t number = quartetToNumber(i, j, k, l);
					numberToQuartet(number, a, b, c, d);
					std::cout << "Quartet (" << i << "," << j << "|" << k << "," << l << ") has number " << number
							<< std::endl;
					std::cout << "But number " << number << " has quartet (" << a << "," << b << "|" << c << "," << d
							<< ")" << std::endl;
				}
			}
		}
	}
}
