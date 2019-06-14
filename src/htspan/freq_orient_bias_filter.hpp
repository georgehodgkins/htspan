#include <gsl/gsl_cdf.h>

#include "math.hpp"
#include "base_orient_bias_filter.hpp"
#include "orient_bias_data.hpp"


#ifndef _HTSPAN_FREQ_ORIENT_BIAS_HPP_
#define _HTSPAN_FREQ_ORIENT_BIAS_HPP_

namespace hts {

using namespace std;

struct freq_orient_bias_filter_f : public base_orient_bias_filter_f {

	freq_orient_bias_filter_f (orient_bias_data &dref)
		: base_orient_bias_filter_f (dref)
	{
	}

	// pass through extra parameters to the base class
	freq_orient_bias_filter_f (orient_bias_data &dref, 
		double lb = -15.0, double ub = 15.0, double eps = 1e-6, size_t max_iter = 100)
		: base_orient_bias_filter_f(dref, lb, ub, eps, max_iter),
	{
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
		double theta_0 = estimate_initial_theta();
		double theta_hat;
		if (known_phi) {
			theta_hat = estimate_theta_given(phi, theta_0);
		} else {
			theta_and_phi theta_phi = estimate_theta_phi(theta_0, phi);
			theta_hat = theta_phi.theta;
			phi = theta_phi.phi;
		}
		double dev = deviance_theta(theta_hat, phi);
		return 1 - gsl_cdf_chisq_P(dev, 1);
	}

};

}// namespace hts

#endif // _HTSPAN_FREQ_ORIENT_BIAS_HPP_