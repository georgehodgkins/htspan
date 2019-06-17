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
#include "functor.hpp"
#include "numeric_integration.hpp"


namespace hts {

using namespace std;


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
	* @template Integrator Integrator object to use, which has an integrate method with signature (double) (double, double, double)
	* @param alpha Value of alpha for the prior beta distribution of phi
	* @param beta Value of beta for the prior beta distribution of phi
	* @return Struct containing evidence values for null and alt models
	*/
	template <typename Integrator>
	evidences model_evidence(double alpha, double beta) {
		// evaluate null model evidence (theta = 0)
		lp_bases_theta_phi_f p_f (*this, alpha, beta, 0);
		Integrator integrator;
		double ev_null = integrator.integrate(p_f, 0, 1);
		// evaluate alternate model evidence (integrating across possible values of theta)
		lp_bases_theta_f t_f (grid_phi, *this, alpha, beta);
		double ev_alt = integrator.integrate(t_f, 0, 1);
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
	template <typename Integrator>
	double operator () (double prior_alt, double alpha, double beta) {
		evidences ev = model_evidence<Integrator>(alpha, beta);
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

	/**
	* Log probability of observed bases given phi and the set theta,
	* plus(+) the log(pdf) of the beta distribution defined by 
	* the set alpha and beta, at the given phi.
	*/
	struct lp_bases_theta_phi_f : public math::numeric_functor {
		// pointer to class containing the lp_bases_given fcn
		bayes_orient_bias_filter_f &bobfilter;
		// alpha for the beta distribution
		double alpha;
		// beta for the beta distribution
		double beta;
		// theta for the log probability of variant function
		double theta;
		// constructor which sets hyperparameters
		lp_bases_theta_phi_f (bayes_orient_bias_filter_f &fi, double a, double b, double t) :
			bobfilter(fi), alpha(a), beta(b), theta(t) {}
		// evaluates the function to be integrated
		double operator() (double phi) {
			return exp(bobfilter.lp_bases_given(theta, phi) +
				log(gsl_ran_beta_pdf(phi, alpha, beta)));
		}
	};


	/**
	* Numerical integration of phi_integrand across a given phi space.
	*/
	struct lp_bases_theta_f : public math::numeric_functor {
		// grid of points which define rectangles for midpoint quadrature
		vector<double> grid_phi;
		// pointer to class containing the lp_bases_given fcn
		bayes_orient_bias_filter_f &bobfilter;
		// alpha for the beta distribution
		double alpha;
		// beta for the beta distribution
		double beta;
		// constructor which sets hyperparameters
		lp_bases_theta_f (vector<double> g, bayes_orient_bias_filter_f &fi, double a, double b) :
			grid_phi(g), bobfilter(fi), alpha(a), beta(b) {}
		// evaluates the function to be integrated (in this case, integral of another function)
		double operator() (double theta) {
			lp_bases_theta_phi_f p_f (bobfilter, alpha, beta, theta);
			math::midpoint integrator;
			return integrator.integrate(p_f, 0, 1);
		}
	};

}; // bayes_orient_bias_filter_f struct

// Integrand functors which take a double x-value and have hyperparameters set at instantiation
// Base class defined in functor.hpp (just an empty virtual operator())

} // namespace hts

#endif // _HTSPAN_BAYES_ORIENT_BIAS_HPP_