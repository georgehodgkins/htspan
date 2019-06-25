#ifndef _HTSPAN_BASE_ORIENT_BIAS_HPP_
#define _HTSPAN_BASE_ORIENT_BIAS_HPP_


#include <cassert>
#include <stdexcept>

#include <string.h>

#include "math.hpp"
#include "orient_bias_data.hpp"
#include "functor.hpp"
#include "brent.hpp"
#include "nucleotide.hpp"

namespace hts {

using namespace std;

// return struct for coordinate ascent estimation
struct theta_and_phi {
	double theta;
	double phi;
};

/**
* This class is the base for the functors which
* implement various damage identification methods.
* It cannot be implemented directly.
*
* This class contains a reference to an underlying data object
* which should be populated before derived classes are instantiated.
* It also contains the methods which implement the damage-adjusted
* variant identfication model, and methods which provide initial 
* estimates for the theta and phi parameters to its children.
*/

struct base_orient_bias_filter_f {
	const char* text_id;// "orientation" set in constructor
	//reference to externally initialized data object
	const orient_bias_data& data;

	//upper and lower bounds for estimation
	const double minimizer_lb;//default -15.0
	const double minimizer_ub;//default 15.0
	//maximum error threshold for estimation (and integration in Bayesian method)
	const double epsabs;//default 1e-6

	base_orient_bias_filter_f(const orient_bias_data &dref,
		double lb = -15.0, double ub = 15.0, double eps = 1e-6, size_t max_iter = 100)
	: text_id("orientation"),
		data(dref),	
		minimizer_lb(lb),
		minimizer_ub(ub),
		epsabs(eps)
	{
	}

	size_t size() const {
		return data.size();
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
			} else if (data.bases[r] == 1) {
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
	double estimate_theta_given(double phi) {
		nlp_bases_given_theta_f f (*this, phi);
		//call common minimization code
		double lmin_theta = math::argmin(f, minimizer_lb, minimizer_ub, epsabs);
		return logistic(lmin_theta);
	}
	
	/**
	* Find value of phi that maximizes the log probability of observed bases,
	* given a value of theta.
	*/
	double estimate_phi_given(double theta) {
		// instantiate functor of objective function to minimize
		nlp_bases_given_phi_f f (*this, theta);
		// call common minimization code
		double lmin_phi = math::argmin(f, minimizer_lb, minimizer_ub, epsabs);
		return logistic(lmin_phi);
	}

	/**
	* Return a human-readable description of what this filter detects.
	* Used for VCF annotation.
	*/
	const char* get_description () const {
		if (data.r1_alt == nuc_T) {
			if (data.r1_ref == nuc_G) {
				return "damage artifact detected by oxoG orientation bias filter";
			}
			if (data.r1_ref == nuc_C) {
				return "damage artifact detected by FFPE orientation bias filter";
			}
		} else {
			return "damage artifact detected by other orientation bias filter";
		}
		// this point should not be reached but it makes the compiler happy
		return NULL;
	}


	struct nlp_bases_given_theta_f : public math::numeric_functor {
		// reference to class containing lp_bases_given
		base_orient_bias_filter_f &filter;
		// fixed phi for theta estimation
		double phi;
		// constructor to set hyperparameters
		nlp_bases_given_theta_f (base_orient_bias_filter_f &fi, double p) :
			filter(fi), phi(p) {}
		// negative log probability of bases given theta
		double operator() (double theta_real) const {
			 return - filter.lp_bases_given(logistic(theta_real), phi);
		}
	};

	struct nlp_bases_given_phi_f : public math::numeric_functor {
		// reference to class containing lp_bases_given
		base_orient_bias_filter_f &filter;
		// fixed theta for phi estimation
		double theta;
		// constructor to set hyperparameters
		nlp_bases_given_phi_f (base_orient_bias_filter_f &fi, double t) :
			filter(fi), theta(t) {}
		// negative log probability of bases given phi
		double operator() (double phi_real) const {
			return - filter.lp_bases_given(theta, logistic(phi_real));
		}
	};
 
};


}  // namespace hts

#endif  // _HTSPAN_BASE_ORIENT_BIAS_HPP_
