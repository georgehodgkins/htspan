#ifndef _HTSPAN_ORIENT_BIAS_QUANT_HPP_
#define _HTSPAN_ORIENT_BIAS_QUANT_HPP_

#include <iostream>
#include <queue>
#include <vector>
#include <cassert>

#include <stdint.h>
#include <string.h>

#include <htslib/hts.h>
#include <htslib/faidx.h>

#include "bam.hpp"
#include "math.hpp"
#include "piler.hpp"
#include "io/faidx_reader.hpp"

namespace hts {

using namespace std;

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

	// reader for BAM sequence data
	piler &p;

	// reader for FASTA sequence data (reference)
	faidx_reader &faidx;

	/**
	 * Initialize class.
	 *
	 * Initialize class with expected ref and alt nucleotides for read 1
	 */
	orient_bias_quant_f(nuc_t _ref, nuc_t _alt)
	: r1_ref(_ref),
		r1_alt(_alt),
		r2_ref(nuc_complement(_ref)),
		r2_alt(nuc_complement(_alt)),
		xc(0), xi(0), nc(0), ni(0)
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
	bool push(bam1_t* b, int32_t pos) {
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

	int next () {
		const bam_pileup1_t *pile = p.next();
		if (pile == NULL) {
			return 1;
		}
		const char* seq = faidx.get(p.tid, p.pos, p.pos+1);
		if (seq == NULL) {
			return 2;
		}
		// check that reference nucleotide is consistent
		char s_ref = char_to_nuc(seq[0]);
		if (s_ref == r1_ref || s_ref == r2_ref) {
			read_count += push(pile, p.size(), p.pos);
		} else {
			return 3;
		}
		return 0;
	}

	virtual size_t push (const bam_pileup1_t *pile, size_t n, int32_t pos) = 0;

	virtual size_t size() const = 0;

	virtual double operator()() const = 0;

	//TODO: add simulation options for testing

};

struct freq_orient_bias_quant_f : base_orient_bias_quant_f {

	// count of read sites
	size_t site_count;

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
			if (push(pile[i].b, pos)) {
				++success;
			}
		}
		if (success > 0) {
			++site_count;
		}
		return success;
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


struct bayes_orient_bias_quant_f : base_orient_bias_quant_f {

	// vector of consistent alt counts, by site
	vector<long int> xc_vec;

	// vector of inconsistent alt counts, by site
	vector<long int> xi_vec;

	// vector of consistent total counts, by site
	vector<long int> nc_vec;

	//vector of inconsistent total counts, by site
	vector<long int> ni_vec;

	size_t push(const bam_pileup1_t* pile, size_t n, int32_t pos) {
		xc = 0;
		xi = 0;
		nc = 0;
		ni = 0;
		size_t success = 0;
		for (size_t i = 0; i < n; ++i) {
			if (push(pile[i].b, pos)) {
				++success;
			}
		}
		xc_vec.push_back(xc);
		xi_vec.push_back(xi);
		nc_vec.push_back(nc);
		ni_vec.push_back(ni);
	}

	size_t size() const {
		return xc_vec.size();
	}
	
	struct theta_hparams_optimizable {

		// parent object containing data
		const bayes_orient_bias_quant_f &parent;

		// vector of parameters updated during optimization
		// [0] = alpha_theta
		// [1] = beta_theta
		vector<double> curr_params;

		theta_hparams_optimizable (bayes_orient_bias_quant_f &P, vector<double> &initial_params) :
			parent(P), curr_params(initial_params) {}

		size_t nobs() const {
			return parent.size();
		}

		size_t nparams() const {
			return 2;
		}

		static double lp_xij_given_hparams (const long int xij, const long int nij, double alpha_theta, double beta_theta) {
			return lbeta(alpha_theta + xij, beta_theta + nij - xij) + lchoose(nij, xij); 
		}

		// AKA -lp_xi_given_hparams
		double operator() (const vector<double> &x) const {
			// xi_vec and ni_vec are class members
			// x[0] is alpha_theta, x[1] is beta_theta
			double sum_lp_xij = 0.0;
			for (size_t j = 0; j < nobs(); ++j) {
				sum_lp_xij += lp_xij_given_hparams(parent.xi_vec[j], parent.ni_vec[j], x[0], x[1]);
			}
			return -(sum_lp_xij - lbeta(x[0], x[1])*nobs());
		}

		void accumulate (vector<double> &current_grad) {
			vector<double> new_grad (2);
			stograd::finite_difference_gradient(*this, curr_params, new_grad);
			stograd::add_to(new_grad, current_grad);
		}

		void update (vector<double> &delta) {
			stograd::subtract_from(delta, curr_params);
		}
	}

	struct phi_hparams_optimizable {

		// parent object containing data
		const bayes_orient_bias_quant_f &parent;

		const double alpha_theta;
		const double beta_theta;

		// vector of parameters updated during optimization
		// [0] = alpha_phi
		// [1] = beta_phi
		vector<double> curr_params;

		phi_hparams_optimizable(const bayes_orient_bias_quant_f &P, vector<double> &initial_params,
				const double a_theta, const double b_theta) :
				parent(P), curr_params(initial_params), alpha_theta(a_theta), beta_theta(b_theta) {}

		size_t nobs() const {
			return parent.size();
		}

		size_t nparams() const {
			return 2;
		}

		static double lp_xcj_given_hparams (const long int xcj, const long int ncj, const double alpha_theta, const double beta_theta,
				const double alpha_phi, const double beta_phi) {
			double *lse_array = new double[xcj+1];
			for (int k = 0; k <= xcj; ++k) {
				lse_array[k] =
					-log(ncj - k + 1) +
					lbeta(xcj - k + alpha_phi, ncj - xcj + beta_phi) -
					lbeta(xcj - k + 1,         ncj - xcj + 1) +
					lbeta(k + alpha_theta, ncj - k + beta_theta) -
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
			for (size_t j = 0; j < nobs(); ++j) {
				sum_lp_xcj += lp_xcj_given_hparams(parent.xc_vec[j], parent.nc_vec[j], alpha_theta, beta_theta, x[0], x[1]);
				sum_log_ncp += log(nc + 1);
			}
			return -(sum_lp_xcj - sum_log_ncp - nobs()*(lbeta(alpha_theta, beta_theta) + lbeta(x[0], x[1])));
		}

		void accumulate (vector<double> &current_grad) {
			vector<double> new_grad (2);
			stograd::finite_difference_gradient(*this, curr_params, new_grad);
			stograd::add_to(new_grad, current_grad);
		}

		void update (vector<double> &delta) {
			stograd::subtract_from(delta, curr_params);
		}
	}
};

}  // namespace hts

#endif  // _HTSPAN_ORIENT_BIAS_QUANT_HPP_
