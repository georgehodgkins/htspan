// This header contains helper functions and test thresholds used in unit testing

#ifndef _HTSPAN_TEST_HPP_
#define _HTSPAN_TEST_HPP_

#include <cmath>

#ifndef TEST_EPS
#define TEST_EPS 1e-3
#endif

#ifndef TEST_RAT
#define TEST_RAT .05
#endif

// Returns the difference ratio of two values, the distance between the values divided by their mean
inline double diff_rat (double a, double b) {
	return 2.0*abs(a-b)/abs(a+b);
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
		return diff_rat(got, std) < TEST_RAT;
	}
}

// true if value is close to or less than the standard
bool test_val_lte (double got, double std) {
	return (got <= std) || test_val(got, std);
}

#endif // _HTSPAN_TEST_HPP