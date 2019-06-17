#ifndef _HTSPAN_OPTIMIZATION_HPP_
#define _HTSPAN_OPTIMIZATION_HPP_

#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_cdf.h>

#include "functor.hpp"

namespace hts {

inline bool point_is_below_endpoints(gsl_function *f, double x, double f_lb, double f_ub) {
	return (f->function(x, f->params) < f_lb && f->function(x, f->params) < f_ub);
}

/**
* This function tries some other guesses if the initial guess is not in the bounds.
* Should be replaced with proper meta-analysis (Bayesian optimization).
*/
double naive_point_picker (gsl_function *f, double x_0, double lb, double ub, size_t max_pow = 4) {
	size_t p = 1;
	double f_lb = f->function(lb, f->params);
	double f_ub = f->function(ub, f->params);
	while (p <= max_pow) {
		double div = pow(2, p);
		double off = (x_0 - lb)/div;
		if (point_is_below_endpoints(f, lb + off, f_lb, f_ub)) {
			return lb + off;
		}
		if (p > 1 && point_is_below_endpoints(f, x_0 - off, f_lb, f_ub)) {
			return x_0 - off;
		}
		off = (ub - x_0)/div;
		if (point_is_below_endpoints(f, x_0 + off, f_lb, f_ub)) {
			return x_0 + off;
		}
		if (p > 1 && point_is_below_endpoints(f, ub - off, f_lb, f_ub)) {
			return ub - off;
		}
		++p;
	}
	return x_0;
}




/**
* Optimization routine, using Brent method implementation in GSL.
*
* CAVEATS:
* If the function value at x_0 is not less than the values at the endpoints,
* the function will return the smaller of the endpoints and will not perform a minimization.
*
* If the guess is not within the range (ub, lb), the minimizer will return the guess unmodified.
*
* This function will only find one minimum in the given bounds, typically the one closest to the lower bound.
* 
* @param f_ptr Pointer to a gsl_function struct initialized by the caller
* @param x_0 An initial guess for the minimizer, within (0, 1)
* @param minimizer_lb Lower bound to minimize on
* @param minimizer_ub Upper bound to minimize on
* @param max_minimizer_iter Maximum number of minimizer iterations
* @param epsabs Threshold for convergence of minimizer
* @return Minimum from the GSL minimizer, if possible (see above)
*/


double minimizer_base (gsl_function *f, double x_0,
		const double minimizer_lb, const double minimizer_ub, const size_t max_minimizer_iter, const double epsabs) {

	// if the function has already been minimized to an endpoint or
	// starts outside the bounds simply return the guess
	if (x_0 <= minimizer_lb || x_0 >= minimizer_ub) {
		return x_0;
	}
	
	// if the function at the initial guess is not lower than
	// the value at both endpoints, minimizer will not work,
	// so we try searching some other areas for a valid point
	double f_lb = f->function(minimizer_lb, f->params);
	double f_ub = f->function(minimizer_ub, f->params);
	double f_x0 = f->function(x_0, f->params);
	if (f_x0 > f_lb || f_x0 > f_ub) {
		return (f_lb < f_ub) ? minimizer_lb : minimizer_ub;
	}
	
	// initialize minimizer
	gsl_min_fminimizer* s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
	gsl_min_fminimizer_set(s, f, x_0, minimizer_lb, minimizer_ub);
	
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
	double xmin = gsl_min_fminimizer_x_minimum(s);
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
	// this point should never be reached but this makes the compiler happy
	return GSL_NAN;
}


/**
* Two wrappers for the above minimization routine,
* which convert numeric_functors to appropriate gsl_function objects.
* 
* Argmin returns the x-value at which the function minimum occurs (+/- eps).
* Argmax returns the x-value at which the fucntion maximum occurs (+/- eps).
*
* See the caveats on minimizer_base above for important details about functionality.
*/
double argmin (numeric_functor &func, double x_0,
		const double minimizer_lb, const double minimizer_ub, const double max_minimizer_iter, const double epsabs) {
	// Calls standard (non-negative) functor conversion method 
	gsl_function f = func.to_gsl_function();
	return minimizer_base(&f, x_0, minimizer_lb, minimizer_ub, max_minimizer_iter, epsabs);
}


double argmax (numeric_functor &func, double x_0,
		const double minimizer_lb, const double minimizer_ub, const double max_minimizer_iter, const double epsabs) {
	// Calls negated functor conversion method (minimizing the negative is maximization)
	gsl_function f = func.to_neg_gsl_function();
	return minimizer_base(&f, x_0, minimizer_lb, minimizer_ub, max_minimizer_iter, epsabs);
}

} // namespace hts

#endif // _HTSPAN_OPTIMIZATION_HPP_