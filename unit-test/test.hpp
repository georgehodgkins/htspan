// This header contains helper functions and test thresholds used in unit testing

#ifndef _HTSPAN_TEST_HPP_
#define _HTSPAN_TEST_HPP_

#include <cmath>

#ifndef TEST_EPS
#define TEST_EPS 1e-6
#endif

#ifndef TEST_RAT
#define TEST_RAT .05
#endif

// Returns the difference ratio of two values, the distance between the values divided by their mean
inline double diff_rat (double a, double b) {
	return 2.0*abs(a-b)/abs(a+b);
}

/**
* Differentiated closeness tests for values of various sizes.
*
* If both values are larger than TEST_EPS, their difference ratio
* must be less than TEST_RAT for them to pass.
*
* Otherwise, if either value is zero, the magnitude of the other value
* must be less than TEST_EPS^2 for them to pass.
*
* If neither value is zero but at least one is very small, they must be within
* an order of magnitude (base 10) to pass.
*/
bool test_val(double got, double std) {
	// true equality
	if (got == std) {
		return true;
	}
	if (abs(got) < TEST_EPS || abs(std) < TEST_EPS) {
		// tests for values extremely close to zero
		if (got == 0.0) {
			return abs(std) < TEST_EPS*TEST_EPS;
		}
		if (std == 0.0) {
			return abs(got) < TEST_EPS*TEST_EPS;
		}
		// one order of magnitude, base 10
		return (abs(got/std) < 10) && (abs(got/std) > .1);
	} else {
		return diff_rat(got, std) < TEST_RAT;
	}
}

// true if value is close to or less than the standard
// Used in minimization (a value less than the standard is better)
bool test_val_lte (double got, double std) {
	return (got <= std) || test_val(got, std);
}

#endif // _HTSPAN_TEST_HPP