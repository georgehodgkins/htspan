#ifndef _HTSPAN_FREQ_ORIENT_BIAS_HPP_
#define _HTSPAN_FREQ_ORIENT_BIAS_HPP_

#include <gsl/gsl_cdf.h>

#include "math.hpp"
#include "base_orient_bias_filter.hpp"
#include "orient_bias_data.hpp"

namespace hts {

using namespace std;

/**
* This class implements a frequentist model for variant identification.
*
* Its operator() returns a p-value for the genuine variant, given an estimate
* for global damage and a specification of whether phi is fixed or should also
* be estimated.
*/

struct freq_orient_bias_filter_f : public base_orient_bias_filter_f {

	// pass through extra parameters to the base class
	freq_orient_bias_filter_f (orient_bias_data &dref, 
		double lb = -15.0, double ub = 15.0, double eps = 1e-6)
		: base_orient_bias_filter_f(dref, lb, ub, eps)
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
	* Estimate theta and phi jointly using coordinate ascent.
	*
	* @param theta initial value of theta
	* @param phi initial value of phi
	* @param eps numerical tolerance
	* @param max_iter maximum number of iterations
	* @return a struct containing the estimated values of theta and phi
	*/
	theta_and_phi estimate_theta_phi(double theta, double phi, int max_iter = 100) {
		double theta_old = theta;
		double phi_old = phi;
		for (int n = 1; n <= max_iter; ++n) {
			theta = estimate_theta_given(phi);
			phi = estimate_phi_given(theta);
			if (abs(theta_old - theta) < epsabs && abs(phi_old - phi) < epsabs) {
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
			theta_hat = estimate_theta_given(phi);
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