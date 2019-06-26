#ifndef _HTSPAN_MATH_HPP_
#define _HTSPAN_MATH_HPP_

#include <cmath>

#include "functor.hpp"
#include "de_integrator.hpp"
#include "gamma.hpp"

namespace hts {

using namespace std;

/**
 * Logit.
 */
double logit(double x) {
	return log(x) - log(1.0 - x);
}

/**
 * Logistic.
 */
double logistic(double x) {
	return 1.0 / (1.0 + exp(-x));
}

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

/**
* This functor represents the function that is integrated
* when calculating the beta function: t^(a-1) * (1-t)^(b-1)
*/
struct beta_kernel_f : math::numeric_functor {
	// hyperparameters
	double alpha;
	double beta;
	// constructor which sets hyperparameters
	beta_kernel_f (double a, double b) :
		alpha(a), beta(b) {}
	double operator() (double t) const {
		return pow(t, alpha - 1.0)*pow(1.0 - t, beta - 1.0);
	}
};

double lbeta (double alpha, double beta) {
	return math::lgamma(alpha) + math::lgamma(beta) - math::lgamma(alpha + beta);
}

/**
* The log probability density of the beta distribution with the
* hyperparameters alpha and beta at x.
*
* PDF = x^(a-1)*(1-x)^(b-1)/B(a, b) where B is the beta fcn:
* definite integral of t^(a-1)*(1-t)^(b-1) on (0, 1)
*
* Then the log PDF is: (a-1)/log(x) + (b-1)/log(1-x) - log(B(a,b))
*
* The log beta function (denominator) is computationally expensive,
* so the hyperparameters are cached and it is only re-evaluated if 
* they change.
*
* @param x Point at which to evaluate the pdf
* @param alpha Alpha hyperparameter of the beta distribution
* @param beta Beta hyperparameter of the beta distribution
* @return The log PDF at x for the beta dist defined by alpha and beta
*/
double log_beta_pdf (double x, double alpha, double beta) {
	// cached parameters and log beta function value
	static double alpha_cached = -1.0;
	static double beta_cached = -1.0;
	static double lbeta_f_cached = -1.0;
	// update log beta function value if parameters changed
	if (alpha != alpha_cached || beta != beta_cached) {
		lbeta_f_cached = lbeta(alpha, beta);
		alpha_cached = alpha;
		beta_cached = beta;
	}
	// evaluate log PDF
	return (alpha_cached - 1.0)*log(x) + (beta_cached - 1.0)*log(1.0 - x) - lbeta_f_cached;
}

}  // namespace hts

#endif  // HTSPAN_MATH_HPP_

