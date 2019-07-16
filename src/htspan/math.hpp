#ifndef _HTSPAN_MATH_HPP_
#define _HTSPAN_MATH_HPP_

#include <cmath>

#include "alglib/specialfunctions.h"

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

// These methods are taken from the OpenCV source at
// modules/hal/include/opencv2/hal/defs.h
bool isinf(double x)
{
    union { uint64_t u; double f; } ieee754;
    ieee754.f = x;
    return ( (unsigned)(ieee754.u >> 32) & 0x7fffffff ) == 0x7ff00000 &&
           ( (unsigned)ieee754.u == 0 );
}

bool isnan(double x)
{
    union { uint64_t u; double f; } ieee754;
    ieee754.f = x;
    return ( (unsigned)(ieee754.u >> 32) & 0x7fffffff ) +
           ( (unsigned)ieee754.u != 0 ) > 0x7ff00000;
}

double fmin2(double x, double y) {
	return (x < y) ? x : y;
}

double fmax2(double x, double y) {
	return (x > y) ? x : y;
}

/**
*  Container-generic element-wise log, modifying
*/
template <typename Container>
Container& log_elements (Container &X) {
	for (typename Container::iterator it = X.begin(); it != X.end(); ++it) {
		*it = log(*it);
	}
	return X;
}

/**
*  Container-generic element-wise log, non-modifying
*/
template <typename Container>
Container log_c (Container X) {
	for (typename Container::iterator it = X.begin(); it != X.end(); ++it) {
		*it = log(*it);
	}
	return X;
}

/**
* Container-generic element-wise exp, modifying
*/
template <typename Container>
Container& exp_elements (Container &X) {
	for (typename Container::iterator it = X.begin(); it != X.end(); ++it) {
		*it = exp(*it);
	}
	return X;
}

/**
* Container-generic element-wise exp, non-modifying
*/
template <typename Container>
Container exp_c (Container X) {
	for (typename Container::iterator it = X.begin(); it != X.end(); ++it) {
		*it = exp(*it);
	}
	return X;
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
* Returns the CDF of the chi-squared distribution
* with n degrees of freedom at the value x.
*/
double chisq_cdf (double x, double n) {
	if (x <= 0.0) {
		return 0;
	}
	return alglib::chisquaredistribution(n, x);
}

/**
* Returns the log of the beta function for 
* parameters alpha and beta.
*/
double lbeta (const double alpha, const double beta) {
	double foo;// return value set by reference that we don't need
	return alglib::lngamma(alpha, foo) + alglib::lngamma(beta, foo) - alglib::lngamma(alpha + beta, foo);
}

// log(n!/(k!(n-k)!))
double lchoose (const int n, const int k) {
	double nfac_log = 0.0;// log(n!)
	double kfac_log = 0.0;// log(k!)
	double nkdfac_log = 0.0;// log((n-k)!)
	int G; // max(k, n-k)
	if (n-k >= k) {
		G = n-k;
		for (int x = 2; x <= k; ++x) {
			kfac_log += log(x);
		}
		nkdfac_log = kfac_log;
		for (int x = k+1; x <= n-k; ++x) {
			nkdfac_log += log(x);
		}
		nfac_log = nkdfac_log;
	} else {
		G = k;
		for (int x = 2; x <= n-k; ++x) {
			nkdfac_log += log(x);
		}
		kfac_log = nkdfac_log;
		for (int x = n-k+1; x <= k; ++x) {
			kfac_log += log(x);
		}
		nfac_log = kfac_log;
	}
	for (int x = G+1; x <= n; ++x) {
		nfac_log += log(x);
	}
	return nfac_log - kfac_log - nkdfac_log;
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

