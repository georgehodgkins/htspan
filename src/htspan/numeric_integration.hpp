#ifndef _HTSPAN_NUMERIC_INTEGRATION_HPP_
#define _HTSPAN_NUMERIC_INTEGRATION_HPP_

// This header file defines various functors and methods for numeric integration

namespace hts {

using namespace std;

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
vector<double> generate_cgrid (double center, double eps = 1e-6, double step = .005, double base = 1.4) {
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
	vector<double> grid (v_cap);
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
* will have to be set externally (here, typically as class members)
*
* @param grid Vector of x-values which define grid.size()-1 rectangles to use for integration
* @param f Integrand object, which has an operator() method that takes and returns double.
*/
double midpoint_integration (vector<double> grid, base_integrand_f &f) {
	// calculate midpoints and grid widths
	vector<double> midpoints (grid.size() - 1);
	vector<double> widths (grid.size() - 1);
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

}// namespace hts

#endif // _HTSPAN_NUMERIC_INTEGRATION_HPP_