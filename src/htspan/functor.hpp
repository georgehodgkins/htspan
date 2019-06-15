#ifndef _HTSPAN_FUNCTOR_HPP_
#define _HTSPAN_FUNCTOR_HPP_

#include <gsl/gsl_math.h>

namespace hts {

// forward declaration of functor type for the following functions
template<typename Real> struct numeric_functor;

/**
* Evaluates a numeric_functor at the given value of x.
*/
template <typename Real>
Real evaluate_numeric_functor (Real x, void* v) {
	numeric_functor<Real> *p = (numeric_functor<Real>*) v;
	return (*p)(x);
}

/**
* Returns the negative of a functor's value at the given x.
*/
template <typename Real>
inline Real neg_evaluate_numeric_functor(Real x, void* v) {
	return - evaluate_numeric_functor(x, v);
}

/**
* This struct acts as a base interface for functors
* passed to numerical methods, which must implement
* an operator() method which takes double and returns double.
* The only other method in derived classes should be a constructor
* which sets integrand hyperparameters.
*
* Note that while this type is templated for compatibility with
* Boost class signatures, it is not designed to work with types other than double,
* and non-double instantiations of the type are not supported.
*
* Compatible with Boost quadrature methods.
*/
template <typename Real>
struct numeric_functor {
	virtual Real operator() (Real x) {
		return 0.0;
	}
	// return gsl_function corresponding to this functor
	gsl_function to_gsl_function () {
		gsl_function f;
		f.function = &evaluate_numeric_functor<Real>;
		f.params = (void*) this;
		return f;
	}
	// return negated gsl_function corresponding to this functor
	gsl_function to_neg_gsl_function () {
		gsl_function f;
		f.function = &neg_evaluate_numeric_functor<Real>;
		f.params = (void*) this;
		return f;
	}

};

} // namespace hts

#endif // _HTSPAN_FUNCTOR_HPP_