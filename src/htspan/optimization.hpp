#ifndef _HTSPAN_OPTIMIZATION_HPP_
#define _HTSPAN_OPTIMIZATION_HPP_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_cdf.h>

#include "numeric_integration.hpp" 

namespace hts {

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
* @param x_0 An initial guess for the minimizer, within (0, 1)
* @return Minimum from the GSL minimizer, if possible (see above)
*/

double minimize_log_function (numeric_functor &func, double x_0) {
	// create function object to pass to minimizer
	gsl_function f;
	f.function = &evaluate_numeric_functor;
	f.params = (void*) &func;

	// transform the guess into unconstrained space
	double lx_0 = logit(x_0);
	
	// if the function has already been minimized to an endpoint or
	// starts outside the bounds simply return the guess
	// generally lx_0 will only be OOB if it's -inf (x_0 == 0)
	if (lx_0 <= minimizer_lb || lx_0 >= minimizer_ub) {
		return x_0;
	}
	
	// if the function at the initial guess is not lower than
	// the value at both endpoints, minimizer will not work
	double f_x0 = f_ptr->function(lx_0, f->params);
	double f_lb = f_ptr->function(minimizer_lb, f->params);
	double f_ub = f_ptr->function(minimizer_ub, f->params);
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
	// this point should never be reached but this makes the compiler happy
	return GSL_NAN;
}

} // namespace hts

#endif // _HTSPAN_OPTIMIZATION_HPP_