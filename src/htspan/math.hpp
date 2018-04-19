#ifndef _HTSPAN_MATH_HPP_
#define _HTSPAN_MATH_HPP_

#include <cmath>

namespace hts {

using namespace std;

/**
 * Phred transformation
 */
double phred(double x) {
	return -10 * log10(x);
}

/**
 * Anti-Phred transformation
 */
double anti_phred(double x) {
	return pow(10, -x/10);
}

/**
 * Calculate log(sum(exp(xs))) of array xs.
 *
 * @param n   size of array
 * @param xs  array
 * @return log(sum(exp(xs)))
 */
double log_sum_exp(size_t n, double xs[]) {
	if (n > 0) {
		// find maximum value
		double max = xs[0];
		for (size_t i = 1; i < n; ++i) {
			if (xs[i] > max) {
				max = xs[i];
			}
		}

		// calculate log sum exp
		double sum = 0;
		for (size_t i = 0; i < n; ++i) {
			sum += exp(xs[i] - max);
		}
		return log(sum) + max;
	} else {
		return 0;
	}
}

}  // namespace hts

#endif  // HTSPAN_MATH_HPP_

