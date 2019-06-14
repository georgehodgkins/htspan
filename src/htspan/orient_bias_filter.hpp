#ifndef _HTSPAN_ORIENT_BIAS_HPP_
#define _HTSPAN_ORIENT_BIAS_HPP_


//TODO: cleanup (remove redundant headers, non specific math fcns)
#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <stdexcept>

#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_cdf.h>

#include "math.hpp"

namespace hts {

using namespace std;

inline double _nlp_bases_given_theta(double theta_real, void* p);
inline double _nlp_bases_given_phi(double phi_real, void* p);

struct theta_and_phi {
	double theta;
	double phi;
};

struct orient_bias_filter_f {
	//reference to externally initialized data object
	const orient_bias_data& data;

	// global damage probability
	double phi_t;
	
	// alternate allele frequency (temporary)
	double theta_t;

	//upper and lower bounds for estimation
	const double minimizer_lb;//default -15.0
	const double minimizer_ub;//default 15.0
	//maximum error threshold for estimation
	const double epsabs;//default .001
	//maximum number of iterations for estimation
	const size_t max_minimizer_iter;//default 100

	orient_bias_filter_f(const orient_bias_data &dref,
		double lb = -15.0, double ub = 15.0, double eps = 1e-6, size_t max_iter = 100)
	: data(dref),
		phi_t(0.0),
		theta_t(0.0),	
		minimizer_lb(lb),
		minimizer_ub(ub),
		epsabs(eps),
		max_minimizer_iter(max_iter)
	{
	}

	size_t size() const {
		return bases.size();
	}

	/**
	 * Log probability of reference base given parameters and read r at locus j.
	 * 
	 * p( b_{j,r} = \tt{R} \mid \theta_j, phi_{j,r}, e_{j,r} )
	 */
	static double lp_ref_base_given(double theta, double phi_r, double e_r) {
		const size_t T = 3;
		double lxs[T] = {
			log(theta)     + log(e_r)-log(3) + log(1 - phi_r),
			log(1 - theta) + log(1 - e_r)    + log(1 - phi_r),
			log(theta)     + log(1 - e_r)    + log(phi_r)
		};
		return log_sum_exp(T, lxs);
	}

	/**
	 * Log probability of alternative base given parameters and read r at locus j.
	 * 
	 * p( b_{j,r} = \tt{A} \mid \theta_j, phi_{j,r}, e_{j,r} )
	 */
	static double lp_alt_base_given(double theta, double phi_r, double e_r) {
		const size_t T = 3;
		double lxs[T] = {
			log(theta)     + log(1 - e_r)    + log(1 - phi_r),
			log(1 - theta) + log(e_r)-log(3) + log(1 - phi_r),
			log(1 - theta) + log(1 - e_r)    + log(phi_r)
		};
		return log_sum_exp(T, lxs);
	}

	/**
	 * Log probability of other base given parameters and read r at locus j.
	 * 
	 * p( b_{j,r} = \tt{O} \mid \theta_j, phi_{j,r}, e_{j,r} )
	 */
	static double lp_other_base_given(double theta, double phi_r, double e_r) {
		const size_t T = 2;
		double lxs[T] = {
			log(e_r)-log(3) + log(1 - phi_r),
			log(1 - e_r)    + log(phi_r)
		};
		return log_sum_exp(T, lxs);
	}
	
	/**
	* Estimate initial value of theta from reads.
	*/
	double estimate_initial_theta() {
		double theta_0 = 1.0;
		for (size_t r = 0; r < data.bases.size(); ++r) {
			if (data.bases[r] == 1) {
				++theta_0;
			}
		}
		theta_0 /= (data.bases.size() + 2.0);
		return theta_0;
	}
	
	/**
	 * Log probability of observed bases given parameters.
	 *
	 */
	double lp_bases_given(double theta, double phi) const {
		double lp = 0;
		// marginalize over all reads.
		size_t R = data.bases.size();
		for (size_t r = 0; r < R; ++r) {
			double phi_r = data.orients[r] ? phi : 0;
			if (data.bases[r] == 0) {
				lp += lp_ref_base_given(theta, phi_r, data.errors[r]);
			} else if (bases[r] == 1) {
				lp += lp_alt_base_given(theta, phi_r, data.errors[r]);
			} else {
				lp += lp_other_base_given(theta, phi_r, data.errors[r]);
			}
		}
		return lp;
	}

	/**
	 * Find value of theta that maximizes the log probability of observed bases,
	 * given a value of phi.
	 */
	double estimate_theta_given(double phi, double theta_0) {
		// define objective function to minimize
		phi_t = phi;
		gsl_function f;
		f.function = &_nlp_bases_given_theta;
		f.params = this;
		//call common minimization code
		return minimize_log_function(&f, theta_0);
	}
	
	/**
	*  Find value of phi that maximizes the log probability of observed bases,
	* given a value of theta.
	*/ 
	double estimate_phi_given(double theta, double phi_0) {
		// define objective function to minimize
		theta_t = theta;
		gsl_function f;
		f.function = &_nlp_bases_given_phi;
		f.params = this;
		//call common minimization code
		return minimize_log_function(&f, phi_0);
	}

	/**
	* Minimization code common to the phi and theta estimation functions, using Brent method.
	*
	* If the function value at x_0 is not less than the values at the endpoints,
	* the function will return the smaller of the endpoints and will not perform a minimization.
	*
	* epsabs, minimizer_ub, minimizer_lb, and max_minimizer_iter
	* are const class members set in constructor
	* 
	* @param f_ptr Pointer to a gsl_function struct initialized by the caller
	* @param x_0 An initial guess for the minimizer
	* @return Minimum from the GSL minimizer, if possible (see above)
	*/

	double minimize_log_function (gsl_function *f_ptr, double x_0) {
		double lx_0 = logit(x_0); // transform the guess into unconstrained space
		// if the function has already been minimized to an endpoint or
		// starts outside the bounds simply return the guess
		// generally lx_0 will only be OOB if it's -inf (x_0 == 0)
		if (lx_0 <= minimizer_lb || lx_0 >= minimizer_ub) {
			return x_0;
		}
		// if the function at the initial guess is not lower than
		// the value at both endpoints, minimizer will not work
		double f_x0 = f_ptr->function(lx_0, this);
		double f_lb = f_ptr->function(minimizer_lb, this);
		double f_ub = f_ptr->function(minimizer_ub, this);
		if (f_x0 > f_lb || f_x0 > f_ub) {
			return (f_lb < f_ub) ? logistic(minimizer_lb) : logistic(minimizer_ub);
		}
		// initialize minimizer
		gsl_min_fminimizer* s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
		gsl_min_fminimizer_set(s, f_ptr, lx_0, minimizer_lb, minimizer_ub);
		// iterate on the minimizer
		int status = GSL_CONTINUE;
		for (size_t iter = 0; iter < max_minimizer_iter; ++iter) {
			status = gsl_min_fminimizer_iterate(s);
			if (status == GSL_FAILURE || status == GSL_EBADFUNC) break;
			// test for success
			double a = gsl_min_fminimizer_x_lower(s);
			double b = gsl_min_fminimizer_x_upper(s);
			status = gsl_min_test_interval(a, b, epsabs, 0.0);
			if (status == GSL_SUCCESS) break;
		}
		// store minimizer result and then deallocate
		double xmin = logistic(gsl_min_fminimizer_x_minimum(s));
		gsl_min_fminimizer_free(s);
		switch (status) {
			case GSL_SUCCESS:
				return xmin;
			case GSL_FAILURE:
			case GSL_EBADFUNC:
				throw runtime_error("Unknown error in minimizer");
			default:
				throw runtime_error("Minimizer did not converge");
		}
		//this point should never be reached but this makes the compiler happy
		return GSL_NAN;
	}

	/**
	* Estimate theta and phi jointly using coordinate ascent.
	*
	* @param theta initial value of theta
	* @param phi initial value of phi
	* @param eps numerical tolerance
	* @param max_iter maximum number of iterations
	* @return a struct containing the estimated values of theta and phi
	*/
	theta_and_phi estimate_theta_phi(double theta, double phi, double eps = 1e-6, int max_iter = 100) {
		double theta_old = theta;
		double phi_old = phi;
		for (int n = 1; n <= max_iter; ++n) {
			theta = estimate_theta_given(phi, theta);
			phi = estimate_phi_given(theta, phi);
			if (abs(theta_old - theta) < eps && abs(phi_old - phi) < eps) {
				break;
			}
			theta_old = theta;
			phi_old = phi;
		}
		theta_and_phi rtn;
		rtn.theta = theta;
		rtn.phi = phi;
		return rtn;
	}

	/**
	 * Calculate deviance of fitted and null models.
	 *
	 * model_1: theta = theta.hat
	 * model_0: theta = 0
	 */
	double deviance_theta(double theta_hat, double phi) {
		return 2 * (lp_bases_given(theta_hat, phi) - lp_bases_given(0.0, phi));
	}

	/**
	 * Variant test with adjustment for orientation bias.
	 *
	 * Test whether variant allele frequence > 0,
	 * with consideration for orietation bias, as specified
	 * by the estimated global damage rate causing  orientation bias.
	 *
	 * Relevant reads should have been tallied using push() already.
	 *
	 * @param phi estimated global damage rate
	 * @param known_phi consider phi known/fixed (true) or variable (false)
	 * @return p-value
	 */
	double operator()(double phi, bool known_phi=true) {
		volatile double theta_0 = estimate_initial_theta();
		volatile double theta_hat;
		if (known_phi) {
			theta_hat = estimate_theta_given(phi, theta_0);
		} else {
			theta_and_phi theta_phi = estimate_theta_phi(theta_0, phi);
			theta_hat = theta_phi.theta;
			phi = theta_phi.phi;
		}
		volatile double dev = deviance_theta(theta_hat, phi);
		return 1 - gsl_cdf_chisq_P(dev, 1);
	}


// theta_real is in unconstrained, real-number space
inline double _nlp_bases_given_theta(double theta_real, void* p) {
	orient_bias_filter_f* x = (orient_bias_filter_f*) p;
	return - x->lp_bases_given(logistic(theta_real), x->phi_t);
}


// phi_real is in unconstrained, real-number space
inline double _nlp_bases_given_phi(double phi_real, void* p) {
	orient_bias_filter_f* x = (orient_bias_filter_f*) p;
	return - x->lp_bases_given(x->theta_t, logistic(phi_real));
}

}  // namespace hts

#endif  // _HTSPAN_ORIENT_BIAS_HPP_
