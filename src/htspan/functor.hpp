#ifndef _HTSPAN_FUNCTOR_HPP_
#define _HTSPAN_FUNCTOR_HPP_

namespace hts {

/**
* This struct acts as a base interface for functors
* passed to numerical methods, which must implement
* an operator() method which takes double and returns double.
* The only other method in derived classes should be a constructor
* which sets integrand hyperparameters.
*
* Compatible with Boost quadrature methods.
*/
struct numeric_functor {
	virtual double operator() (double x) {}
}

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
inline double neg_evaluate_numeric_functor(double x, void* v) {
	return - evaluate_numeric_functor(x, v);
}

} // namespace hts

#endif // _HTSPAN_FUNCTOR_HPP_