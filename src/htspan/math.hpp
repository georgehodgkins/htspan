#ifndef _HTSPAN_MATH_HPP_
#define _HTSPAN_MATH_HPP_

#include <cmath>
#include <algorithm>

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
	double foo;
	return alglib::lngamma(alpha, foo) + alglib::lngamma(beta, foo) - alglib::lngamma(alpha + beta, foo);
}

// log(n!/(k!(n-k)!))
double lchoose (const int n, const int k) {
	double nfac_log = 0.0;//log(n!)
	double kfac_log = 0.0;//log(k!)
	double nkdfac_log = 0.0;//log((n-k)!)
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
* Generates a batch of binomially distributed random numbers
* by generating a u~(0, 1) and naively finding
* the largest k such that u < Binom(k, n, p).
*
* @param n Number of trials
* @param p Probability of success in each trial
* @param n_vars Number of variates to generate
* @return a vector of variates <= n
*/
void rbinom_batch(int n, double p, size_t n_vars, vector<int> &K) {
	K.reserve(n_vars);
	// generate n_vars number of u~(0,1)
	vector<double> U (n_vars);
	alglib::hqrndstate rng;
	alglig::hqrndrandomize(rng);
	for (size_t i = 0; i < n_vars; ++i) {
		U.push_back(alglib::hqrnduniformr(rng));
	}
	// sort in descending order
	sort(U.rbegin(), U.rend());
	int k = 0;
	// evaluate the distribution at each value of k, advancing to the next u-value each time it increases past the previous one
	double q = 1 - p;
	double B = pow(q, n);
	while (!U.empty()) {
		do {
			++k;
			B *= (n - k + 1)/k * p/q;
		} while (B < U.back());
		U.pop_back();
		K.push_back(k);
	}
	// shuffle the return vector, since it will be sorted after the above process
	for (size_t i = 0; i < K.size(); ++i) {
		size_t si = alglib::hqrnduniformi(rng, K.size());
		swap(K[i], K[si]);
	}
	return K;
}

/**
* Single random binomial variate for 
* n trials with probability of success p.
*/
// TODO: make this function not suck
int rbinom (int n, double p) {
	static alglib::hqrndstate rng;
	alglib::hqrndrandomize(rng);
	double u = alglib::hqrnduniformr(rng);
	double B = pow(1 - p, n);
	int x = 0;
	while (u < B) {
		++x;
		B *= (n - x + 1)/x * p/q;
	}
	return x;
}

/**
* Single beta-distributed variate for shape parameters
* a and b.
*
* Algorithm is from Cheng, R.C.H., CACM April 1978.
*/
// TODO: add optimizations
double rbeta (double a, double b) {
	static double alpha = 0.0;
	static double beta = 0.0;
	static double gamma = 0.0;
	static double a_cached = -1.0;
	static double b_cached = -1.0;
	if (a != a_cached || b != b_cached) {
		a_cached = a;
		b_cached = b;
		alpha_cached = a + b;
		if (a <= 1 || b <= 1) {
			beta = (1/a < 1/b) ? 1/a : 1/b;
		} else {
			beta = sqrt((alpha - 2)/(2*a*b-alpha));
		}
		gamma = a + 1/beta;
	}
	static alglib::hqrndstate rng;
	alglib::hqrndrandomize(rng);
	double u1, u2, W, V;
	do {
		u1 = alglib::hqrnduniformr(rng);
		u2 = alglib::hqrnduniformr(rng);
		V = beta*log(u1/(1-u1));
		W = a*exp(V);
	} while (alpha*log(alpha/(b+W)) + gamma*V - log(4) < log(u1*u1*u2));
	return W/(b + W);
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

