#ifndef _HTSPAN_ORIENT_BIAS_QUANT_HPP_
#define _HTSPAN_ORIENT_BIAS_QUANT_HPP_

#include <iostream>
#include <vector>
#include <cassert>
#include <cfloat>

#include <stdint.h>
#include <string.h>

#include <htslib/hts.h>
#include <htslib/faidx.h>

#include <alglib/ap.h>
#include <alglib/alglibmisc.h>
#include <alglib/specialfunctions.h>

#include <stograd/src/stograd/stograd.hpp>

#include "bam.hpp"
#include "math.hpp"
#include "r_rand.hpp"
#include "piler.hpp"
#include "io/faidx_reader.hpp"

namespace hts {

// The stepper to use for stochastic gradient optimization in the Bayes model
typedef stograd::stepper::adam<double> bayes_stepper_t;

using namespace std;

struct hparams {
	double alpha_theta;
	double beta_theta;
	double alpha_phi;
	double beta_phi;
};

struct base_orient_bias_quant_f {

	/// reference and alternative nucleotides to consider
	nuc_t r1_ref, r1_alt, r2_ref, r2_alt;

	// oxoG-consistent alt count
	long int xc;
	
	// oxoG-inconsistent alt count
	long int xi;

	// oxoG-consistent total count
	long int nc;
	
	// oxoG-inconsistent total count
	long int ni;

	// count of total reads pushed
	int read_count;

	/**
	 * Initialize class.
	 *
	 * @param _ref Read 1 reference nucleotide for SNV of interest
	 * @param _alt Read 2 reference nucleotide for SNV of interest
	 * @param _p Initialized hts::piler containing sequence data of interest
	 * @param _f Initialized hts::faidx_reader containing reference sequence
	 */
	base_orient_bias_quant_f(nuc_t _ref, nuc_t _alt)
	: r1_ref(_ref),
		r1_alt(_alt),
		r2_ref(nuc_complement(_ref)),
		r2_alt(nuc_complement(_alt)),
		xc(0), xi(0), nc(0), ni(0), read_count(0)
	{
	}

	size_t n_reads() const {
		return read_count;
	}

	/**
	 * Push a read to accumulate statistics.
	 *
	 * @param b       BAM record of a read
	 * @param pos     reference position
	 */
	bool push(const bam1_t* b, int32_t pos) {
		// only analyze nucleotides A, C, G, T (no indels)
		nuc_t qnuc = query_nucleotide(b, pos);
		if (!nuc_is_canonical(qnuc)) return false;
		// query aligning against the reverse strand is reverse-complemented,
		// so we need to reverse-complement again to get the original nucleotide
		nuc_t onuc;
		if (bam_is_rev(b)) {
			onuc = nuc_complement(qnuc);
		} else {
			onuc = qnuc;
		}

		// accmulate statistics based on read1 vs. read2 and original nucleotide
		if (bam_is_read1(b)) {
			if (onuc == r1_ref) {
				++nc;
			} else if (onuc == r1_alt) {
				++nc;
				++xc;
			} else if (onuc == r2_ref) {
				++ni;
			} else if (onuc == r2_alt) {
				++ni;
				++xi;
			}
			// other nucleotides are ignored
		// double-check that the read is second read, in case the flag is malformed
		} else if (bam_is_read2(b)) {
			if (onuc == r2_ref) {
				++nc;
			} else if (onuc == r2_alt) {
				++nc;
				++xc;
			} else if (onuc == r1_ref) {
				++ni;
			} else if (onuc == r1_alt) {
				++ni;
				++xi;
			}
			// other nucleotides are ignored
		} else {
			// read is neither read1 or read2
			// BAM file is likely not from paired end sequencing
			return false;
		}

		return true;
	}

	static void simulate_orient_bias_read_counts (const vector<size_t> &ns_vec,
			const vector<double> &theta_vec, const vector<double> &phi_vec,
			vector<long int> &xc_vec, vector<long int> &xi_vec, vector<long int> &nc_vec, vector<long int> &ni_vec) {
		// setup return vectors
		xc_vec.clear();
		xc_vec.reserve(ns_vec.size());
		xi_vec.clear();
		xi_vec.reserve(ns_vec.size());
		nc_vec.clear();
		nc_vec.reserve(ns_vec.size());
		ni_vec.clear();
		ni_vec.reserve(ns_vec.size());
		for (size_t j = 0; j < ns_vec.size(); ++j) {
			ni_vec.push_back(ns_vec[j]/2);
			nc_vec.push_back(ns_vec[j] - ni_vec[j]);
			xi_vec.push_back(r_rand::rbinom(ni_vec[j], theta_vec[j]));
			long int xr = r_rand::rbinom(nc_vec[j], theta_vec[j]);
			xc_vec.push_back(xr + r_rand::rbinom(nc_vec[j] - xr, phi_vec[j])); // xd ~ Binom(nc, phi) and xc = xr + xd
		}
	}

	/**
	* Accumulate statistics for a set of reads at the given locus.
	*
	* @param pile Pointer to HTSlib pileup object populated with reads of interest
	* @param n Number of reads to process
	* @param pos Reference position of the given pileup
	* @return Number of reads pushed.
	*/
	virtual size_t push (const bam_pileup1_t *pile, size_t n, int32_t pos) = 0;

	/**
	* Return the number of processed loci
	* (n_reads returns total processed reads).
	*/
	virtual size_t size() const = 0;

};

struct freq_orient_bias_quant_f : public base_orient_bias_quant_f {

	// count of read sites
	size_t site_count;

	freq_orient_bias_quant_f(nuc_t _ref, nuc_t _alt) : 
			base_orient_bias_quant_f(_ref, _alt),
			site_count(0)
		{
		}

	size_t size () const {
		return site_count;
	}

	/**
	 * Push reads to accumulate statistics.
	 *
	 * @param pile bam pileup object
	 * @param n    number of reads in pileup object
	 * @param pos  target reference position
	 * @return number of successfully processed reads
	 */
	size_t push(const bam_pileup1_t* pile, size_t n, int32_t pos) {
		size_t success = 0;
		for (size_t i = 0; i < n; ++i) {
			if (base_orient_bias_quant_f::push(pile[i].b, pos)) {
				++success;
			}
		}
		if (success > 0) {
			++site_count;
		}
		return success;
	}

	void simulate(size_t N, double theta, double phi) {
		vector<size_t> ns_vec (1, N);
		vector<double> theta_vec (1, theta);
		vector<double> phi_vec (1, phi);
		vector<long int> xc_vec;
		vector<long int> xi_vec;
		vector<long int> nc_vec;
		vector<long int> ni_vec;
		simulate_orient_bias_read_counts(ns_vec, theta_vec, phi_vec, xc_vec, xi_vec, nc_vec, ni_vec);
		xc = xc_vec.front();
		xi = xi_vec.front();
		nc = nc_vec.front();
		ni = ni_vec.front();
	}

	/**
	 * Return estimate of global DNA damage
	 *
	 * Reads should have been processed by calling `push`.
	 *
	 * @return estimate of phi
	 */
	double operator()() const {
		double theta_hat = double(xi)/ni;
		double phi = (double(xc)/nc - theta_hat) / (1 - theta_hat);
		if (phi < 0.0) {
			return 0.0;
		}
		return phi;
	}
};


struct bayes_orient_bias_quant_f : public base_orient_bias_quant_f {

	// vector of consistent alt counts, by site
	vector<long int> xc_vec;

	// vector of inconsistent alt counts, by site
	vector<long int> xi_vec;

	// vector of consistent total counts, by site
	vector<long int> nc_vec;

	//vector of inconsistent total counts, by site
	vector<long int> ni_vec;

	bayes_orient_bias_quant_f(nuc_t _ref, nuc_t _alt) :
			base_orient_bias_quant_f(_ref, _alt),
			xc_vec(0), xi_vec(0), nc_vec(0), ni_vec(0)
		{
		}

	size_t push(const bam_pileup1_t* pile, size_t n, int32_t pos) {
		xc = 0;
		xi = 0;
		nc = 0;
		ni = 0;
		size_t success = 0;
		for (size_t i = 0; i < n; ++i) {
			if (base_orient_bias_quant_f::push(pile[i].b, pos)) {
				++success;
			}
		}
		xc_vec.push_back(xc);
		xi_vec.push_back(xi);
		nc_vec.push_back(nc);
		ni_vec.push_back(ni);
		return success;
	}

	size_t size() const {
		return xc_vec.size();
	}

	void simulate(const vector<size_t> &ns_vec, const double alpha_theta, const double beta_theta,
			const double alpha_phi, const double beta_phi) {
		vector<double> theta_vec (ns_vec.size());
		vector<double> phi_vec(ns_vec.size());
		for (size_t j = 0; j < ns_vec.size(); ++j) {
			theta_vec[j] = r_rand::rbeta(alpha_theta, beta_theta);
		}
		for (size_t j = 0; j < ns_vec.size(); ++j) {
			phi_vec[j] = r_rand::rbeta(alpha_phi, beta_phi);
		}
		simulate_orient_bias_read_counts(ns_vec, theta_vec, phi_vec, xc_vec, xi_vec, nc_vec, ni_vec);
	}

	hparams operator()(size_t bsize, size_t nepochs, double learning_rate, double eps,
			double alpha0_theta, double beta0_theta, double alpha0_phi, double beta0_phi) {
		vector<double> theta_init(2);
		theta_init[0] = alpha0_theta;
		theta_init[1] = beta0_theta;
		theta_hparams_optimizable theta_opt (*this, theta_init);
		bayes_stepper_t stepper (learning_rate);
		stograd::optimize(theta_opt, stepper, bsize, nepochs, eps);
		vector<double> phi_init(2);
		phi_init[0] = alpha0_phi;
		phi_init[1] = beta0_phi;
		phi_hparams_optimizable phi_opt (*this, phi_init, theta_opt.alpha(), theta_opt.beta());
		stograd::optimize(phi_opt, stepper, bsize, nepochs, eps);
		hparams rtn;
		rtn.alpha_theta = theta_opt.alpha();
		rtn.beta_theta = theta_opt.beta();
		rtn.alpha_phi = phi_opt.alpha();
		rtn.beta_phi = phi_opt.beta();
		return rtn;
	}

	double eval_theta_objective_func (hparams x, int maxJ = -1) {
		vector<double> X (2);
		X[0] = x.alpha_theta;
		X[1] = x.beta_theta;
		theta_hparams_optimizable F (*this, X);
		if (maxJ == -1) {
			F.J = size()-1;
		} else {
			F.J = maxJ;
		}
		return F(X);
	}

	double eval_theta_grad_dalpha (hparams x, int maxJ = -1) {
		vector<double> X (2);
		X[0] = x.alpha_theta;
		X[1] = x.beta_theta;
		theta_hparams_optimizable F(*this, X);
		if (maxJ == -1) {
			F.J = size()-1;
		} else {
			F.J = maxJ;
		}
		return F.dlp_xi_given_hparams_dalpha(x.alpha_theta, x.beta_theta);
	}

	double eval_theta_grad_dbeta (hparams x, int maxJ = -1) {
		vector<double> X (2);
		X[0] = x.alpha_theta;
		X[1] = x.beta_theta;
		theta_hparams_optimizable F(*this, X);
		if (maxJ == -1) {
			F.J = size()-1;
		} else {
			F.J = maxJ;
		}
		return F.dlp_xi_given_hparams_dbeta(x.alpha_theta, x.beta_theta);
	}

	double eval_phi_objective_func (hparams x, int maxJ = -1) {
		vector<double> X (2);
		X[0] = x.alpha_phi;
		X[1] = x.beta_phi;
		phi_hparams_optimizable F(*this, X, x.alpha_theta, x.beta_theta);
		if (maxJ == -1) {
			F.J = size()-1;
		} else {
			F.J = maxJ;
		}
		return F(X);
	}

	void set_data (vector<long int> &new_xc, vector<long int> &new_xi, vector<long int> &new_nc, vector<long int> &new_ni) {
		xc_vec.clear();
		xc_vec = new_xc;
		xi_vec.clear();
		xi_vec = new_xi;
		nc_vec.clear();
		nc_vec = new_nc;
		ni_vec.clear();
		ni_vec = new_ni;
	}

	struct theta_hparams_optimizable {

		// parent object containing data
		const bayes_orient_bias_quant_f &parent;

		// index of record (in parent object) currently being analyzed
		size_t J;

		// vector of parameters being optimized
		// [0] = alpha_theta
		// [1] = beta_theta
		vector<double> curr_params;

		double alpha() const {
			return curr_params[0];
		}

		double beta() const {
			return curr_params[1];
		}

		theta_hparams_optimizable (bayes_orient_bias_quant_f &P, vector<double> &initial_params) :
			parent(P), J(0), curr_params(initial_params) {}

		size_t nobs() const {
			return parent.xi_vec.size();
		}

		size_t cobs() const {
			return J+1;
		}

		size_t nparams() const {
			return 2;
		}

		size_t next () {
			if (J == parent.size()-1) {
				J = 0;
			} else {
				++J;
			}
			return J;
		}

		static double lp_xij_given_hparams (const long int xij, const long int nij, const double alpha_theta, const double beta_theta) {
			return lbeta(alpha_theta + xij, beta_theta + nij - xij) + lchoose(nij, xij); 
		}

		static double dlp_xij_given_hparams_dalpha (const long int xij, const long int nij, const double alpha_theta, const double beta_theta) {
			return alglib::psi(alpha_theta + xij) - alglib::psi(alpha_theta + beta_theta + nij);
		}

		static double dlp_xij_given_hparams_dbeta (const long int xij, const long int nij, const double alpha_theta, const double beta_theta) {
			return alglib::psi(beta_theta + nij - xij) - alglib::psi(alpha_theta + beta_theta + nij);
		}

		double dlp_xi_given_hparams_dalpha (const double alpha_theta, const double beta_theta) {
			static double alpha_cached = -1.0;
			static double psi_alpha_cached = 0.0;
			if (alpha_theta != alpha_cached) {
				psi_alpha_cached = alglib::psi(alpha_theta);
				alpha_cached = alpha_theta;
			}
			// xi_vec and ni_vec are class members
			double sum_dlp_xij_dalpha = 0.0;
			for (size_t j = 0; j < cobs(); ++j) {
				sum_dlp_xij_dalpha += dlp_xij_given_hparams_dalpha(parent.xi_vec[j], parent.ni_vec[j], alpha_theta, beta_theta);
			}
			return sum_dlp_xij_dalpha + cobs()*(alglib::psi(alpha_theta + beta_theta) - psi_alpha_cached);
		}

		double dlp_xi_given_hparams_dbeta (const double alpha_theta, const double beta_theta) {
			static double beta_cached = -1.0;
			static double psi_beta_cached = 0.0;
			if (beta_theta != beta_cached) {
				psi_beta_cached = alglib::psi(beta_theta);
				beta_cached = beta_theta;
			}
			double sum_dlp_xij_dbeta = 0.0;
			for (size_t j = 0; j < cobs(); ++j) {
				sum_dlp_xij_dbeta += dlp_xij_given_hparams_dbeta(parent.xi_vec[j], parent.ni_vec[j], alpha_theta, beta_theta);
			}
			return sum_dlp_xij_dbeta + cobs()*(alglib::psi(alpha_theta + beta_theta) - psi_beta_cached);
		}

		// AKA -lp_xi_given_hparams
		double operator() (const vector<double> &x) const {
			// xi_vec and ni_vec are class members
			// x[0] is alpha_theta, x[1] is beta_theta
			double sum_lp_xij = 0.0;
			for (size_t j = 0; j < cobs(); ++j) {
				sum_lp_xij += lp_xij_given_hparams(parent.xi_vec[j], parent.ni_vec[j], x[0], x[1]);
			}
			return -(sum_lp_xij - lbeta(x[0], x[1])*cobs());
		}

		void accumulate (vector<double> &current_grad) {
			vector<double> new_grad (2);
			new_grad[0] = dlp_xi_given_hparams_dalpha(curr_params[0], curr_params[1]);
			new_grad[1] = dlp_xi_given_hparams_dbeta(curr_params[0], curr_params[1]);
			stograd::add_to(new_grad, current_grad);
			next();
		}

		void update (vector<double> &delta) {
			stograd::subtract_from(delta, curr_params);
		}
	}; // struct theta_hparams_optimizable

	struct phi_hparams_optimizable {

		// parent object containing data
		const bayes_orient_bias_quant_f &parent;

		// fixed parameters for objective function (previously estimated)
		const double alpha_theta;
		const double beta_theta;

		// vector of parameters being optimized
		// [0] = alpha_phi
		// [1] = beta_phi
		vector<double> curr_params;

		double alpha() const {
			return curr_params[0];
		}

		double beta() const {
			return curr_params[1];
		}

		// index of record (in parent object) currently being analyzed
		size_t J;

		phi_hparams_optimizable(const bayes_orient_bias_quant_f &P, vector<double> &initial_params,
				const double a_theta, const double b_theta) :
				parent(P), alpha_theta(a_theta), beta_theta(b_theta), curr_params(initial_params), J(0) {}

		size_t nobs() const {
			return parent.xi_vec.size();
		}

		size_t cobs() const {
			return J+1;
		}

		size_t nparams() const {
			return 2;
		}

		size_t next() {
			if (J == parent.size() - 1) {
				J = 0;
			} else {
				++J;
			}
			return J;
		}

		static double lp_xcj_given_hparams (const long int xcj, const long int ncj, const double alpha_theta, const double beta_theta,
				const double alpha_phi, const double beta_phi) {
			double *lse_array = new double[xcj+1];
			for (int k = 0; k <= xcj; ++k) {
				lse_array[k] =
					-log(ncj - k + 1) +
					lbeta(alpha_phi + xcj - k, beta_phi + ncj - xcj) -
					lbeta(xcj - k + 1,         ncj - xcj + 1) +
					lbeta(alpha_theta + k, beta_theta + ncj - k) -
					lbeta(k + 1, ncj - k + 1);
			}
			double rtn = log_sum_exp(xcj+1, lse_array);
			delete[] lse_array;
			return rtn;
		}

		// AKA -lp_xc_given_params 
		double operator() (const vector<double> &x) const {
			double sum_lp_xcj = 0.0;
			double sum_log_ncp = 0.0;
			for (size_t j = 0; j < cobs(); ++j) {
				sum_lp_xcj += lp_xcj_given_hparams(parent.xc_vec[j], parent.nc_vec[j], alpha_theta, beta_theta, x[0], x[1]);
				sum_log_ncp += log(parent.nc_vec[j] + 1);
			}
			return -(sum_lp_xcj - sum_log_ncp - cobs()*(lbeta(alpha_theta, beta_theta) + lbeta(x[0], x[1])));
		}

		void accumulate (vector<double> &current_grad) {
			vector<double> new_grad;
			stograd::finite_difference_gradient(*this, curr_params, new_grad);
			stograd::add_to(new_grad, current_grad);
			next();
		}

		void update (vector<double> &delta) {
			stograd::subtract_from(delta, curr_params);
		}
	};// struct phi_hparams_optimizable
};// struct bayes_orient_bias_quant_f

}  // namespace hts

#endif  // _HTSPAN_ORIENT_BIAS_QUANT_HPP_
