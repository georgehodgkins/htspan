#ifndef _HTSPAN_OPTIMIZATION_HPP_
#define _HTSPAN_OPTIMIZATION_HPP_

#include <cmath>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_cdf.h>

#include "functor.hpp"
#include "brent.hpp"

namespace math {

/**
* Initial rough minimization optionally performed before a guess is passed to the GSL minimizer.
* This is useful because the GSL minimizer only accepts guesses with y-values lower than the endpoints,
* and converges much faster for good guesses.
*
* Tail-recursive algorithm divides the space between the given bounds into div # of sections, evaluates
* each internal section boundary (/not/ including the bounds themselves), and picks the minimum evaluation 
* as the new guess. The new bounds are the next and previous section boundaries (possibly including one of the
* outer bounds). These new bounds and guess are then passed to the next iteration, for iter+1 number of iterations.
*
* @param f Function of interest
* @param x_0 Initial guess
* @param lb Lower bound of space to be sectioned
* @param ub Upper bound of space to be sectioned
* @param div Number of sections to divide search space into
* @param iter Remaining number of iterations (after the current one)
* @param margin Guesses closer than this value to the bounds will be discarded
*/
double improve_guess (gsl_function *f, double x_0, const double lb, const double ub,
		const size_t div, const size_t iter, const double margin = 1e-6) {
	double inc = abs(lb-ub)/div;
	if (x_0 < lb + margin) {
		x_0 = lb + margin;
	} else if (x_0 > ub - margin) {
		x_0 = ub - margin;
	}
	double min_g = x_0;
	double min_fval = f->function(x_0, f->params);
	for (size_t n = 1; n < div-1; ++n) {
		double fval = f->function(lb + n*inc + margin, f->params);
		if (fval < min_fval) {
			min_fval = fval;
			min_g = lb + n*inc + margin;
		}
	}
	if (iter == 0) {
		return min_g;
	} else {
		return improve_guess(f, min_g, min_g - inc + margin, min_g + inc - margin, div, iter-1);
	}
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
* @param improve Whether to run the initial rough minimization [true]
* @return Minimum from the GSL minimizer, if possible (see above)
*/

double minimizer_base (gsl_function *f, double x_0,
		const double minimizer_lb, const double minimizer_ub, const size_t max_minimizer_iter,
		const double epsabs, bool improve = true) {

	// if the function has already been minimized to an endpoint or
	// starts outside the bounds simply return the guess
	if (x_0 <= minimizer_lb || x_0 >= minimizer_ub) {
		return x_0;
	}
	// a hacky initial minimization algorithm that's not too expensive
	//if (improve) {
	//	x_0 = improve_guess(f, x_0, minimizer_lb, minimizer_ub, 4, 3);
	//}
	double f_lb = f->function(minimizer_lb, f->params);
	double f_ub = f->function(minimizer_ub, f->params);
	double f_x0 = f->function(x_0, f->params);
	// the epsilon here compensates for cases where the guesser places x_0 extremely close to a bound
	if (f_x0 + epsabs > f_lb || f_x0 + epsabs > f_ub) {
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
			throw std::runtime_error("Unknown error in minimizer");
		default:
			throw std::runtime_error("Minimizer did not converge");
	}
	// this point should never be reached but this makes the compiler happy
	return GSL_NAN;
}


/**
* Two wrappers for the above minimization routine.
* 
* Argmin returns the x-value at which the function minimum occurs (+/- eps).
* Argmax returns the x-value at which the fucntion maximum occurs (+/- eps).
* 
*/
double argmin (numeric_functor &f, double x_0,
		const double minimizer_lb, const double minimizer_ub, const double max_minimizer_iter, const double epsabs) {
	double x_min;
	// x_min is passed by reference ans set to the argmin
	local_min(minimizer_lb, minimizer_ub, epsabs, f, x_min);
	return x_min;
}

double argmax (numeric_functor &f, double x_0,
		const double minimizer_lb, const double minimizer_ub, const double max_minimizer_iter, const double epsabs) {
	// Constructs a negated functor around the passed functor
	negated_functor nf (f);
	double x_min;
	// x_min is passed by reference and set to the argmin of the negative, aka the argmax
	local_min(minimizer_lb, minimizer_ub, epsabs, nf, x_min);
	return x_min;
}

} // namespace hts

#endif // _HTSPAN_OPTIMIZATION_HPP_