#include <vector>
#include <cmath>

#include <gsl/gsl_cdf.h>

#include "orient_bias_filter.hpp"
#include "nucleotide.hpp"


namespace hts {

using namespace std;

// forward declarations of static integrand functions
static double phi_integrand (double phi, void *v);
static double theta_integrand (double theta, void *v);

/**
* Return struct for the model_evidence method in bayes_orient_bias_filter_f,
* containing evidence for null and alternative models.
*/
struct evidence_rtn {
	double null;
	double alt;
};

struct bayes_orient_bias_filter_f : public orient_bias_filter_f {

	// parameters for the beta distribution
	double alpha_phi, beta_phi;

	// grid to integrate phi on (stored here to simplify parameter passing)
	vector<double> grid_phi;

// TODO: make a base class which the Bayesian and non-Bayesian classes both inherit from

	bayes_orient_bias_filter_f (nuc_t _ref, nuc_t _alt, size_t n)
		: orient_bias_filter_f (_ref, _alt, n),
		alpha_phi(0.0), beta_phi(0.0),
		grid_phi(0) {}

	/**
	* Generate a vector of exponentially spaced points between zero and one
	* to numerically integrate using the midpoint method. Points will be
	* distributed, from the center, in steps of step*base^k, where k starts
	* at zero and increases by 1 for each step.
	*
	* The set's smallest and largest members will be eps and 1-eps, respectively.
	* It is also guaranteed to contain the center; with the possible exception
	* of the endpoints, all points in the set will be symmetrical around this center.
	*
	* @param center The x-value which the grid will be symmetrical around
	* @param eps Defines the endpoints of the grid, [eps, 1-eps]
	* @param step Coefficient of the exponential term in step generation (see above)
	* @param base Base of the exponential term in step generation (see above)
	* @return A std::vector containing the grid points
	*/
	std::vector<double> generate_cgrid (double center, double eps = 1e-6, double step = .005, double base = 1.4) {
		// calculate number of grid points to allocate
		// log_b is log base b
		// c+s*b^k < 1 - eps ---> k < log_b ((1-eps-c)/2), k is integer
		size_t right_points = (size_t) floor(log((1-eps-center)/step)/log(base));
		// c-s*b^k > eps ---> k < log_b((c-eps)/s), " "
		size_t left_points = (size_t) floor(log((center - eps)/step)/log(base));
		// three extra points: center and two endpoints
		size_t v_cap = right_points + left_points + 3;
		// allocate vector
		std::vector<double> grid (v_cap);
		// place fixed points
		grid[0] = eps;
		grid[v_cap - 1] = 1 - eps;
		// first point is eps, so last generated point on left is at index left_points
		grid[left_points + 1] = center;
		// generate left points
		for (size_t k = 0; k < left_points; ++k) {
			grid[left_points - k] = step * pow(base, k);
		}
		// generate right points
		for (size_t k = 0; k < right_points; ++k) {
			grid[left_points + 2 + k] = step * pow(base, k);
		}
		return grid;
	}

	/**
	* Custom one-dimensional midpoint integration method which takes a
	* vector of points which define the rectangles to use.
	*
	* Note that any other parameters required for evaluation  of f besides the x-value
	* will have to be set externally (here, typically as class members)
	*
	* @param grid Vector of x-values which define grid.size()-1 rectangles to use for integration
	* @param f Pointer to function to numerically integrate, which takes a double and returns double
	*/
	double midpoint_integration (std::vector<double> grid, double (*f) (double, void*), void *v) {
		// calculate midpoints and grid widths
		std::vector<double> midpoints (grid.size() - 1);
		std::vector<double> widths (grid.size() - 1);
		for (size_t i = 0; i < grid.size() - 1; ++i) {
			widths[i] = grid[i+1] - grid[i];
			midpoints[i] = grid[i] + widths[i]/2;
		}
		// evaluate function at each midpoint, multiply by width, and add to previous total
		double result = 0;
		for (size_t i = 0; i < midpoints.size(); ++i) {
			double evl = f(midpoints[i], v);
			double area = evl * widths[i];
			result += area;
		}
		return result;
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
	evidence_rtn model_evidence(double alpha, double beta) {
		// these class members are accessed by the phi_integrand function
		alpha_phi = alpha;
		beta_phi = beta;
		// generate centered grids to integrate on
		double phi_0 = estimate_phi_given(0, 0.5); // change fixed params?
		// note that grid_phi is a class member
		grid_phi = generate_cgrid(phi_0);// using default params
		double theta_0 = estimate_theta_given(phi_0, 0.5);
		vector<double> grid_theta = generate_cgrid(theta_0);
		// evaluate null model evidence (theta = 0)
		// a single midpoint integration of phi_integrand at theta=0
		double ev_null = theta_integrand(0, (void*) this);
		// evaluate alternate model evidence (integrating across possible values of theta)
		double ev_alt = midpoint_integration(grid_theta, theta_integrand, (void*) this);
		evidence_rtn rtn;
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

}; // functor struct

/**
* Log probability of observed bases given phi
* and an externally set theta_t,
* plus the pdf of the prior beta distribution of phi
* at the given value of phi, using externally set 
* alpha and beta parameters.
*/
static double phi_integrand (double phi, void *v) {
	bayes_orient_bias_filter_f *p = (bayes_orient_bias_filter_f*) v;
	// alpha_phi, beta_phi, and theta_t are class members set externally
	return exp( p->lp_bases_given(p->theta_t, phi) +
		log(gsl_ran_beta_pdf(phi, p->alpha_phi, p->beta_phi)));
}

/**
* Numerical integration of phi_integrand using
* an externally genrated and set grid_phi.
*/
static double theta_integrand (double theta, void *v) {
	bayes_orient_bias_filter_f *p = (bayes_orient_bias_filter_f*) v;
	// grid_phi is a class member set externally
	p->theta_t = theta;
	return p->midpoint_integration(p->grid_phi, phi_integrand, v);
}


} // namespace hts


	