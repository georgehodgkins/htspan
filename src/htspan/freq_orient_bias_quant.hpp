#ifndef _HTSPAN_FREQ_ORIENT_BIAS_QUANT_HPP_
#define _HTSPAN_FREQ_ORIENT_BIAS_QUANT_HPP_

#include "base_orient_bias_quant.hpp"

namespace hts {

struct freq_orient_bias_quant_f : public base_orient_bias_quant_f {

	// previous values of observed variables, used for debug
	long int old_xi, old_xc, old_ni, old_nc;

	// count of read loci (not reads)
	size_t site_count;

	// epsilon for convergence
	const double eps;

	// window size for convergence
	const size_t check_window_size = 10000;

	freq_orient_bias_quant_f(nuc_t _ref, nuc_t _alt, size_t _mx, double _e) : 
			base_orient_bias_quant_f(_ref, _alt, _mx),
			old_xi(0), old_xc(0), old_ni(0), old_nc(0),
			site_count(0), eps(_e)
		{
		}

	/*
	* Returns the count of read loci (not reads).
	* n_reads() returns total reads instead.
	*/
	size_t size () const {
		return site_count;
	}

	void copy_data (const vector<long int> &xc_vec, const vector<long int> &xi_vec,
			const vector<long int> &nc_vec, const vector<long int> &ni_vec) {
		xc = 0;
		xi = 0;
		nc = 0;
		ni = 0;
		for (size_t j = 0; j < xc_vec.size(); ++j) {
			xc += xc_vec[j];
			xi += xi_vec[j];
			nc += nc_vec[j];
			ni += ni_vec[j];
		}
	}

	/**
	 * Process a vector containing
	 * the reads at one locus and update the observed variables accordingly.
	 *
	 * @param pile Already populated vector of pointers to BAM records
	 * @param pos  locus reference position
	 * @return number of successfully processed reads
	 */
	size_t push(const vector<bam1_t*> &pile, int32_t pos) {
		size_t success = 0;

		// copy old values of observed vars
		old_xi = xi;
		old_xc = xc;
		old_ni = ni;
		old_nc = nc;

		// accumulate statistics from new reads
		for (size_t i = 0; i < pile.size(); ++i) {
			if (base_orient_bias_quant_f::push(pile[i], pos)) {
				++success;
			}
		}
		if (success > 0) {
			++site_count;
		}
		read_count += success;
		return success;
	}

	/*
	* Simulate appropriately distributed observed variables, for testing.
	*
	* @param N Total number of reads to simulate
	* @param theta Global alternate allele rate
	* @param phi Global damage rate
	*/
	void simulate(size_t N, double theta, double phi, int seed = 0) {
		simulate_orient_bias_read_counts(N, theta, phi, xc, xi, nc, ni, seed);
	}

	double theta_hat() const {
		return double(xi)/ni;
	}

	/**
	 * Return estimate of global DNA damage (phi_hat).
	 *
	 * Reads should have been processed by calling `push`.
	 *
	 * @return estimate of phi
	 */
	double operator()() const {
		double phi_hat = (double(xc)/nc - theta_hat()) / (1 - theta_hat());
		if (phi_hat < 0.0) {
			return 0.0;
		}
		return phi_hat;
	}

	// return observed variables at last pushed site, for debug
	long int xij () const {return xi - old_xi;}
	long int xcj () const {return xc - old_xc;}
	long int nij () const {return ni - old_ni;}
	long int ncj () const {return nc - old_nc;}

	/*
	* Whether more reads should be pushed.
	* 
	* Returns true if the estimate converges (checked every 10K reads) or
	* the maximum number of reads is reached.
	*
	* Specifically, convergence is checked against the mean of the
	* previous estimates, sampling every 10K reads; possible at 20K reads at the earliest.
	*
	* Also note that if the estimate is less than the set eps, it is not considered convergent
	* (not a big issue; eps should be very small).
	*/
	bool done () const {
		// sum of previously sampled estimates
		static double window_curr = -1;
		static double window_prev = -1;
		static double curr_sum = 0;
		// number of checks performed
		static size_t prev_rdc = 0;
		static size_t n_checks = 1;

		curr_sum += operator()();
		// check for convergence
		bool converged = false;
		// check_window_size is a class constant
		if (read_count > check_window_size*(n_checks + 1)) {
			window_prev = window_curr;
			window_curr = curr_sum / (read_count - prev_rdc);
			prev_rdc = read_count;
			converged = abs(window_prev - window_curr) < eps;
		}
		return converged || base_orient_bias_quant_f::done();
	}

};

} // namespace hts

#endif // _HTSPAN_FREQ_ORIENT_BIAS_QUANT_HPP_
