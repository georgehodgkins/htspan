#ifndef _HTSPAN_BAYES_ORIENT_BIAS_HPP_
#define _HTSPAN_BAYES_ORIENT_BIAS_HPP_

#include <vector>
#include <cmath>
#include <stdexcept>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>

#include "base_orient_bias_filter.hpp"
#include "orient_bias_data.hpp"
#include "nucleotide.hpp"
#include "numeric_integration.hpp"


namespace hts {

using namespace std;

// forward declarations of static integrand functions
static double phi_integrand (double phi, void *v);
static double theta_integrand (double theta, void *v);

/**
* Return struct for the model_evidence method in bayes_orient_bias_filter_f,
* containing evidence for null and alternative models.
*/
struct evidences {
	double null;
	double alt;
};

/**
* This class implements a Bayesian model for variant identification.
*
* Its operator() returns the log posterior probability of a genuine variant
* given the prior variant probability and alpha and beta values for the beta distribution.
*/

struct bayes_orient_bias_filter_f : public base_orient_bias_filter_f {
	//TODO: streamline math methods to use the same calling paradigm

	// parameters for the beta distribution
	double alpha_phi, beta_phi;

	// grid to integrate phi on (stored here to simplify parameter passing)
	vector<double> grid_phi;

	bayes_orient_bias_filter_f (orient_bias_data &dref)
		: base_orient_bias_filter_f (dref),
			alpha_phi(0.0),
			beta_phi(0.0),
			grid_phi(0)
	{
	}

	// pass through extra parameters to the base class
	bayes_orient_bias_filter_f (orient_bias_data &dref, 
		double lb = -15.0, double ub = 15.0, double eps = 1e-6, size_t max_iter = 100)
		: base_orient_bias_filter_f(dref, lb, ub, eps, max_iter),
			alpha_phi(0.0),
			beta_phi(0.0),
			grid_phi(0)
	{
	}

	/**
	* Calculate evidence for the null (theta=0) and alternative
	* (theta > 0) models for orientation bias test.
	*
	* Note that the return type for this method is
	* defined above the main class in this file.
	*
	* @param alpha Value of alpha for the prior beta distribution of phi
	* @param beta Value of beta for the prior beta distribution of phi
	* @return Struct containing evidence values for null and alt models
	*/
	evidences model_evidence(double alpha, double beta) {
		// these class members are accessed by the phi_integrand function
		alpha_phi = alpha;
		beta_phi = beta;
		// generate centered grids to integrate on
		double phi_hat = estimate_phi_given(0, 0.5);
		// note that grid_phi is a class member
		grid_phi = generate_cgrid(phi_hat);// using default params
		double theta_hat = estimate_theta_given(phi_hat, 0.5);
		vector<double> grid_theta = generate_cgrid(theta_hat);
		// evaluate null model evidence (theta = 0)
		// a single midpoint integration of phi_integrand at theta=0
		double ev_null = theta_integrand(0, (void*) this);
		// evaluate alternate model evidence (integrating across possible values of theta)
		double ev_alt = midpoint_integration(grid_theta, theta_integrand, (void*) this);
		evidences rtn;
		rtn.null = ev_null;
		rtn.alt = ev_alt;
		return rtn;
	}

	/**
	* Bayesian model for orientation bias identification.
	* 
	* Computes the posterior probability of the alternative model
	* (theta > 0) against that of the null model (theta = 0).
	*
	* Reads should already have been tallied using push().
	*
	* @param prior_alt Prior probability of the alternative model
	* @param alpha Value of alpha for the prior beta distribution of phi
	* @param beta Value of beta for the prior beta distribution of phi
	* @return Posterior probability of the alternative model
	*/
	double operator () (double prior_alt, double alpha, double beta) {
		evidence_rtn ev = model_evidence(alpha, beta);
		// evaluate the log posterior probability of the alternate model
		double lev_null = log(ev.null);
		double lev_alt = log(ev.alt);

		double lprior_null = log(1-prior_alt);
		double lprior_alt = log(prior_alt);

		double lxs[] = {
			lev_null + lprior_null,
			lev_alt + lprior_alt
		};
		double lse = log_sum_exp(2, lxs);
		double lposterior = lev_alt + lprior_alt - lse;
		return lposterior;
	}

}; // bayes_orient_bias_filter_f struct

// Static integrand functions which take a double x-value and a void pointer to an object containing other parameters
// As used here, the object will be the above class itself

/**
* Log probability of observed bases given phi
* and an externally set theta_t,
* plus the pdf of the prior beta distribution of phi
* at the given value of phi.
*/
static double phi_integrand (double phi, void *v) {
	cgrid_orient_bias_filter_f *p = (cgrid_orient_bias_filter_f*) v;
	// alpha_phi, beta_phi, and theta_t are class members set externally
	return exp( p->lp_bases_given(p->theta_t, phi) +
		log(gsl_ran_beta_pdf(phi, p->alpha_phi, p->beta_phi)));
}

/**
* Numerical integration of phi_integrand across a given phi space.
*/
static double theta_integrand (double theta, void *v) {
	cgrid_orient_bias_filter_f *p = (cgrid_orient_bias_filter_f*) v;
	// grid_phi is a class member set externally
	p->theta_t = theta;
	return p->midpoint_integration(p->grid_phi, phi_integrand, v);
}


} // namespace hts

#endif // _HTSPAN_BAYES_ORIENT_BIAS_HPP_