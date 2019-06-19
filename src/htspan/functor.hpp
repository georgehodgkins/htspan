#ifndef _HTSPAN_FUNCTOR_HPP_
#define _HTSPAN_FUNCTOR_HPP_

#include <gsl/gsl_math.h>

namespace math {

// forward declaration of evaluation functions called in struct methods
//
// Note that the reverse order (declaring struct before functions) does not work
// because the methods of a forward-declared class are not known
double evaluate_numeric_functor (double x, void* v);
double neg_evaluate_numeric_functor(double x, void* v);

/**
* This struct acts as a base interface for functors
* passed to numerical methods, which must implement
* an operator() method which takes a single double and returns double.
* The only other method in derived classes should be a constructor
* which sets integrand hyperparameters.
*
*/
struct numeric_functor {
	// all children of this class must implement an operator()
	// const so it can work with external integrator code
	virtual double operator() (double x) const {
		return 0.0;
	}
	// return gsl_function corresponding to this functor
	gsl_function to_gsl_function () {
		gsl_function f;
		f.function = &evaluate_numeric_functor;
		f.params = (void*) this;
		return f;
	}
	// return negated gsl_function corresponding to this functor
	gsl_function to_neg_gsl_function () {
		gsl_function f;
		f.function = &neg_evaluate_numeric_functor;
		f.params = (void*) this;
		return f;
	}
};

// constructs a functor which returns the negative of the functor passed to it as a parameter
struct negated_functor : public numeric_functor {
	numeric_functor &f;

	negated_functor (numeric_functor &func) : f(func) {}

	double operator() (double x) {
		return - f(x);
	}
};

/**
* Evaluates a numeric_functor at the given value of x.
*/
double evaluate_numeric_functor (double x, void* v) {
	numeric_functor *p = (numeric_functor*) v;
	return (*p)(x);
}

/**
* Returns the negative of a functor's value at the given x.
*/
double neg_evaluate_numeric_functor(double x, void* v) {
	return - evaluate_numeric_functor(x, v);
}

} // namespace hts

#endif // _HTSPAN_FUNCTOR_HPP_