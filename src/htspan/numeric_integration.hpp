#ifndef _HTSPAN_NUMERIC_INTEGRATION_HPP_
#define _HTSPAN_NUMERIC_INTEGRATION_HPP_

// This header file defines various functors and methods for numeric integration

#include <vector>
#include <numeric_limits>

#include "functor.hpp"

namespace hts {

using namespace std;

template <typename Real>
struct midpoint<Real> {

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
	static vector<Real> generate_cgrid (Real center, Real eps = 1e-6, Real step = .005,
			Real base = 1.4, Real lb = 0.0, Real ub = 1.0) {
		// calculate number of grid points to allocate
		// log_b is log base b
		// c+s*b^k < 1 - eps ---> k < log_b ((1-eps-c)/2), k is integer
		size_t right_points = (size_t) ceil(log((1-eps-center)/step)/log(base));
		// c-s*b^k > eps ---> k < log_b((c-eps)/s), " "
		size_t left_points = (size_t) ceil(log((center - eps)/step)/log(base));
		if (center < eps) {
			left_points = 0;
		} else if (center > 1 - eps) {
			right_points = 0;
		}
		// three extra points: center and two endpoints
		size_t v_cap = right_points + left_points + 3;
		// allocate vector
		vector<Real> grid (v_cap);
		// place fixed points
		grid[0] = eps;
		grid[v_cap - 1] = 1 - eps;
		// first point is eps, so last generated point on left is at index left_points
		grid[left_points + 1] = center;
		// generate left points
		for (size_t k = 0; k < left_points; ++k) {
			grid[left_points - k] = center - step * pow(base, k);
		}
		// generate right points
		for (size_t k = 0; k < right_points; ++k) {
			grid[left_points + 2 + k] = center + step * pow(base, k);
		}
		return grid;
	}

	/**
	* One-dimensional midpoint integration method which takes a
	* vector of points which define the rectangles to use.
	*
	* Note that any other parameters required for evaluation  of f besides the x-value
	* are set in the constructor and should not vary during integration
	*
	* @param grid Vector of x-values which define grid.size()-1 rectangles to use for integration
	* @param f numeric_functor object, which has an operator() method that takes and returns template type Real.
	*/
	static Real midpoint_integration (vector<Real> grid, numeric_functor<Real> &f) {
		// calculate midpoints and grid widths
		vector<Real> midpoints (grid.size() - 1);
		vector<Real> widths (grid.size() - 1);
		for (size_t i = 0; i < grid.size() - 1; ++i) {
			widths[i] = grid[i+1] - grid[i];
			midpoints[i] = grid[i] + widths[i]/2;
		}
		// evaluate function at each midpoint, multiply by width, and add to previous total
		Real result = 0;
		for (size_t i = 0; i < midpoints.size(); ++i) {
			Real evl = f(midpoints[i]);
			Real area = evl * widths[i];
			result += area;
		}
		return result;
	}

	/**
	* Autobounding integration method written to have a signature
	* compatible with Boost integrators.
	* 
	* Note that while the method is templated for compatibility with Boost,
	* it is not designed to work with functors other than the numeric_functor type
	* defined in functor.hpp.
	*
	* @param f Functor to integrate, taking a single argument of type Real and returning Real
	* @param a Lower bound of integration
	* @param b Upper bound of integration
	* @param step Step value for integration
	* @return The estimate of the definite integral of F from a to b
	*/
	template <typename F>
	Real integrate (const F f, Real a, Real b, Real step) {
		// passing in the midpoint as a guess is not a good way to do this,
		// since it'll probably kick out the endpoint half the time
		// adding Bayesian optimization would be much better
		Real center = argmax<Real>(f, (a-b)/2, a, b, 50, numeric_limits<Real>::epsilon);
		// generate grid to integrate on
		vector<Real> cgrid = generate_cgrid<Real>(center, numeric_limits<Real>::epsilon, step, 1.4, a, b);
		// do the integration
		return midpoint_integration(cgrid, f);
	}

}// namespace hts

#endif // _HTSPAN_NUMERIC_INTEGRATION_HPP_