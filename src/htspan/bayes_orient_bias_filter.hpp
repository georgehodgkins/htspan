#include <deque>
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
	* Generate a set of exponentially spaced points between zero and one
	* to numerically integrate using the midpoint method. Points will be
	* distributed, from the center, in steps of step*base^x, where x starts
	* at zero and increases by 1 for each step.
	*
	* The set's smallest and largest members will be eps and 1-eps, respectively.
	* It is also guaranteed to contain the center; with the possible exception
	* of the endpoints, all points in the set will be symmetrical around this center.
	*
	* @param center The point which the grid will be symmetrical around
	* @param eps Defines the endpoints of the grid, [eps, 1-eps]
	* @param step Coefficient of the exponential term in step generation (see above)
	* @param base Base of the exponential term in step generation (see above)
	* @return A std::deque containing the grid points
	*/
	std::deque<double> generate_grid (double center, double eps = 1e-6, double step = .005, double base = 1.4) {
		// initialize a deque containing only the center value
		std::deque<double> grid_deck(1, center);
		double exp = 0.0;
		double off = step*pow(base, exp);
		// successively add values to each end of the deque,
		// symmetrically offset from the center
		while (center - off > eps && center + off < 1 - eps) {
			grid_deck.push_back(center + off);
			grid_deck.push_front(center - off);
			++exp;
			off = step*pow(base, exp);
		}
		// ensure that both ends are completed (only one loop will execute)
		while (center - off > eps) {
			grid_deck.push_front(center - off);
			++exp;
			off = step*pow(base, exp);
		}
		while (center + off < 1 - eps) {
			grid_deck.push_back(center + off);
			++exp;
			off = step*pow(base, exp);
		}
		// Add eps and 1-eps as endpoints
		grid_deck.push_front(eps);
		grid_deck.push_back(1 - eps);
		return grid_deck;
	}

	