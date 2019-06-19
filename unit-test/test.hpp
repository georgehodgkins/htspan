// This header contains helper functions and test thresholds used in unit testing

#ifndef _HTSPAN_TEST_HPP_
#define _HTSPAN_TEST_HPP_

#include <cmath>

#define TEST_EPS 1e-3
#define TEST_RAT .5

// Returns the difference ratio of two values, the distance between the values divided by their mean
// not used at the moment
inline double diff_rat (double a, double b) {
	return 2.0*abs(a-b)/(a+b);
}

// Differentiated tests for values smaller than the test eps
// True if the test passes
bool test_val(double got, double std) {
	if (abs(got) < TEST_EPS || abs(std) < TEST_EPS) {
		// tests for values extremely close to zero
		if (got == 0.0) {
			return abs(std) < TEST_EPS*TEST_EPS;
		}
		if (std == 0.0) {
			return abs(got) < TEST_EPS*TEST_EPS;
		}
		return (got/std < 10) && (got/std > .1);
	} else {
		return abs(got-std) < TEST_EPS;
	}
}

#endif // _HTSPAN_TEST_HPP