// This header contains helper functions and test thresholds used in unit testing

#ifndef _HTSPAN_TEST_HPP_
#define _HTSPAN_TEST_HPP_

#define TEST_EPS 1e-4
#define TEST_RAT .5

// Returns the difference ratio of two values, the distance between the values divided by their mean
inline double diff_rat (double a, double b) {
	return 2.0*abs(a-b)/(a+b);
}

// Differentiated tests for values smaller than the test eps
// True if the test passes
bool test_val(double got, double std) {
	if (got < TEST_EPS || std < TEST_EPS) {
		// tests for values extremely close to zero
		if (got == 0.0) {
			return std < TEST_EPS*TEST_EPS;
		}
		if (std == 0.0) {
			return got < TEST_EPS*TEST_EPS;
		}
		return diff_rat(got, std) < TEST_RAT;
	} else {
		return abs(got-std) < TEST_EPS;
	}
}

#endif // _HTSPAN_TEST_HPP