#ifndef _HTSPAN_FUNCTOR_HPP_
#define _HTSPAN_FUNCTOR_HPP_


namespace math {

// TODO: fix unusual behavior if operator() is pure virtual

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
		return 13.337;
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

} // namespace hts

#endif // _HTSPAN_FUNCTOR_HPP_