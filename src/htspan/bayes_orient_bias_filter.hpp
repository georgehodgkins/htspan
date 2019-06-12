#include <vector>
#include <cmath>

#include "orient_bias_filter.hpp"
#include "nucleotide.hpp"

// UNDER CONSTRUCTION

namespace hts {

struct bayes_orient_bias_filter_f : public orient_bias_filter_f {

// TODO: make a base class which the Bayesian and non-Bayesian both inherit from
	bayes_orient_bias_filter_f (nuc_t _ref, nuc_t _alt, size_t n)
		: orient_bias_filter_f (_ref, _alt, n) {}

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
	std::vector<double> generate_grid (double center, double eps = 1e-6, double step = .005, double base = 1.4) {
		// calculate number of grid points to allocate
		// log_b is log base b
		// c+s*b^k < 1 - eps ---> k < log_b ((1-eps-c)/2) 
		size_t right_points = (size_t) floor(log((1-eps-center)/step)/log(base));
		// c-s*b^k > eps ---> k < log_b((c-eps)/s)
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
		for (size_t k = 0; k <= left_points; ++k) {
			grid[left_points - k] = step * pow(base, k);
		}
		// generate right points
		for (size_t k = 0; k <= right_points; ++k) {
			grid[left_points + 2 + k] = step * pow(base, k);
		}
		return grid;
	}

	/**
	* Custom midpoint integration method which takes a
	* vector of points which define the rectangles to use.
	*
	* @param grid Vector of points which define grid.size()-1 rectangles to use for integration
	* @param f Pointer to function to numerically integrate, which takes double and returns double
	*/
	double midpoint_integration (std::vector<double> grid, double (*f) (double)) {
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
			double evl = f(midpoints[i]);
			double area = evl * widths[i];
			result += area;
		}
		return result;
	}
	