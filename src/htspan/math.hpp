#ifndef _HTSPAN_MATH_HPP_
#define _HTSPAN_MATH_HPP_

#include <cmath>

#include "functor.hpp"
#include "de_integrator.hpp"

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
* This functor is the function that is integrated
* when calculating the beta function: t^(a-1) * (1-t)^(b-1)
*/
struct Beta_integrand_f : math::numeric_functor {
	// hyperparameters
	double alpha;
	double beta;
	// constructor which sets hyperparameters
	Beta_fcn_integrand (double a, double b) :
		alpha(a), beta(b) {}
	double operator() (double t) const {
		return pow(t, alpha - 1.0)*pow(1.0 - t, beta - 1.0);
	}
};

/**
* The probability density of the beta distribution with the
* hyperparameters alpha and beta at x.
*
* PDF = x^(a-1)*(1-x)^(b-1)/B(a, b) where B is the beta fcn
*
* The beta function (denominator) is computationally expensive,
* so the hyperparameters are stored and it is only re-evaluated if 
* they change.
*
* @param x Point at which to evaluate the pdf
* @param alpha Alpha hyperparameter of the beta distribution
* @param beta Beta hyperparameter of the beta distribution
* @return The PDF at x for the beta dist defined by alpha and beta
*/
double beta_pdf (double x, double alpha, double beta) {
	// static parameter and Beta fcn evaluation
	static double alpha_stored = -1.0;
	static double beta_stored = -1.0;
	static double B_stored = -1.0;
	// update Beta fcn (denominator) if parameters changed
	if (alpha != alpha_stored || beta != beta_stored) {
		Beta_integrand_f f (alpha, beta);
		math::tanh_sinh<math::numeric_functor> integrator;
		B_stored = integrator.integrate(f, 0, 1, 1e-6);
		alpha_stored = alpha;
		beta_stored = beta;
	}
	// evaluate PDF
	return pow(x, alpha - 1.0)*pow(1.0 - x, beta - 1.0) / B_stored;
}

}  // namespace hts

#endif  // HTSPAN_MATH_HPP_

