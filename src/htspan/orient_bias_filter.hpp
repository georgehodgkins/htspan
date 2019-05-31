#ifndef _HTSPAN_ORIENT_BIAS_HPP_
#define _HTSPAN_ORIENT_BIAS_HPP_

#include <iostream>
#include <queue>
#include <vector>
#include <cassert>

#include <stdint.h>
#include <string.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_cdf.h>

#include "fetcher.hpp"
#include "math.hpp"

namespace hts {

using namespace std;

inline double _nlp_bases_given(double theta_real, void* p);


struct orient_bias_filter_f {

	/// reference base (0), alternative base (1), or other base (2)
	vector<int8_t> bases;
	/// base error probability
	vector<double> errors;
	/// whether variant has an orientation consistent with artifact damage
	/// (false if base is not variant)
	vector<bool> orients;

	// global damage probability
	double phi;

	/// reference and alternative nucleotides to consider
	nuc_t r1_ref, r1_alt, r2_ref, r2_alt;

	/// observed nucleotide of the read
	vector<char> cnucs;
	/// base quality score at the query position
	vector<int> quals;
	/// forward or reverse reference strand to which the read aligned
	vector<char> strands;
	/// first read or second read of a pair
	vector<int> members;

	orient_bias_filter_f(nuc_t _ref, nuc_t _alt, size_t n)
	: r1_ref(_ref),
		r1_alt(_alt),
		r2_ref(nuc_complement(_ref)),
		r2_alt(nuc_complement(_alt)),
		phi(0)
	{
		bases.reserve(n);
		errors.reserve(n);
		orients.reserve(n);
	}

	size_t size() const {
		return bases.size();
	}

	/**
	 * Push a read to accumulate statistics.
	 *
	 * @param b       BAM record of a read
	 * @param pos     reference position
	 * @param nt_ref  reference nucleotide
	 * @param nt_alt  alternative nucleotide
	 */
	bool push(bam1_t* b, int32_t pos, nuc_t nt_ref, nuc_t nt_alt) {
		// only analyze nucleotides A, C, G, T (no indels)
		nuc_t qnuc = query_nucleotide(b, pos);
		if (!nuc_is_canonical(qnuc)) return false;

		bool success = true;

		// determine if query nucleotide matches reference or alternative
		// NB  query aligning against the reverse strand is reverse-complemented
		//     s.t. its sequence is in the same direction as the reference;
		//     the same is true for the reported nt_ref and nt_alt
		if (qnuc == nt_ref) {
			bases.push_back(0);
		} else if (qnuc == nt_alt) {
			bases.push_back(1);
		} else {
			bases.push_back(2);
		}

		// call orientation (damage-consistency) based on the original nucleotide of the read
		nuc_t onuc;
		if (bam_is_rev(b)) {
			onuc = nuc_complement(qnuc);
		} else {
			onuc = qnuc;
		}
		if (bam_is_read1(b)) {
			if (onuc == r1_ref || onuc == r1_alt) {
				// e.g. G or T on read 1
				orients.push_back(true);
			} else {
				orients.push_back(false);
			}
		// double-check that the read is second read, in case the flag is malformed
		} else if (bam_is_read2(b)) {
			if (onuc == r2_ref || onuc == r2_alt) {
				// e.g. C or A on read 2
				orients.push_back(true);
			} else {
				orients.push_back(false);
			}
		} else {
			// flag is malformed
			orients.push_back(false);
			success = false;
		}

		errors.push_back(anti_phred((double)query_quality(b, pos)));

		// populate informational fields
		
		quals.push_back(query_quality(b, pos));

		if (bam_is_rev(b)) {
			strands.push_back('-');
		} else {
			strands.push_back('+');
		}

		if (bam_is_read1(b)) {
			members.push_back(1);
		} else if (bam_is_read2(b)) {
			members.push_back(2);
		} else {
			members.push_back(0);
		}

		switch (qnuc) {
		case nuc_A:
			cnucs.push_back('A');
			break;
		case nuc_C:
			cnucs.push_back('C');
			break;
		case nuc_G:
			cnucs.push_back('G');
			break;
		case nuc_T:
			cnucs.push_back('T');
			break;
		default:
			// this should never happen
			cnucs.push_back('N');
			break;
		}

		if (!success) {
			// rollback everything
			bases.pop_back();
			errors.pop_back();
			orients.pop_back();
			cnucs.pop_back();
			quals.pop_back();
			strands.pop_back();
			members.pop_back();
		}

		return success;
	}

	/**
	 * Push reads to accumulate statistics.
	 *
	 * @param bs   pile of reads
	 * @param pos  target reference position
	 * @param nt_ref  reference nucleotide
	 * @param nt_alt  alternative nucleotide
	 * @return number of successfully processed reads
	 */
	size_t push(vector<bam1_t*> bs, int32_t pos, nuc_t nt_ref, nuc_t nt_alt) {
		size_t success = 0;
		for (size_t i = 0; i < bs.size(); ++i) {
			if (push(bs[i], pos, nt_ref, nt_alt)) {
				++success;
			}
		}
		return success;
	}

	/**
	 * Log probability of reference base given parameters and read r at locus j.
	 * 
	 * p( b_{j,r} = \tt{R} \mid \theta_j, phi_{j,r}, e_{j,r} )
	 */
	static double lp_ref_base_given(double theta, double phi_r, double e_r) {
		const size_t T = 3;
		double lxs[T] = {
			log(theta)     + log(e_r)-log(3) + log(1 - phi_r),
			log(1 - theta) + log(1 - e_r)    + log(1 - phi_r),
			log(theta)     + log(1 - e_r)    + log(phi_r)
		};
		return log_sum_exp(T, lxs);
	}

	/**
	 * Log probability of alternative base given parameters and read r at locus j.
	 * 
	 * p( b_{j,r} = \tt{A} \mid \theta_j, phi_{j,r}, e_{j,r} )
	 */
	static double lp_alt_base_given(double theta, double phi_r, double e_r) {
		const size_t T = 3;
		double lxs[T] = {
			log(theta)     + log(1 - e_r)    + log(1 - phi_r),
			log(1 - theta) + log(e_r)-log(3) + log(1 - phi_r),
			log(1 - theta) + log(1 - e_r)    + log(phi_r)
		};
		return log_sum_exp(T, lxs);
	}

	/**
	 * Log probability of other base given parameters and read r at locus j.
	 * 
	 * p( b_{j,r} = \tt{O} \mid \theta_j, phi_{j,r}, e_{j,r} )
	 */
	static double lp_other_base_given(double theta, double phi_r, double e_r) {
		const size_t T = 2;
		double lxs[T] = {
			log(e_r)-log(3) + log(1 - phi_r),
			log(1 - e_r)    + log(phi_r)
		};
		return log_sum_exp(T, lxs);
	}

	/**
	 * Log probability of observed bases given parameters.
	 *
	 */
	double lp_bases_given(double theta, double phi) const {
		double lp = 0;
		// marginalize over all reads.
		size_t R = bases.size();
		for (size_t r = 0; r < R; ++r) {
			double phi_r = orients[r] ? phi : 0;
			if (bases[r] == 0) {
				lp += lp_ref_base_given(theta, phi_r, errors[r]);
			} else if (bases[r] == 1) {
				lp += lp_alt_base_given(theta, phi_r, errors[r]);
			} else {
				lp += lp_other_base_given(theta, phi_r, errors[r]);
			}
		}
		return lp;
	}

	/**
	 * Find value of theta that maximizes the log probability of observed bases,
	 * given a value of phi.
	 */
	double estimate_theta_given(double phi) {
		const double epsabs = 0.001;
		const size_t max_iter = 100;

		// define objective function to minimize
		this->phi = phi;
		gsl_function f;
		f.function = &_nlp_bases_given;
		f.params = this;

		// initial estimate
		double theta_0 = 1.0;
		for (size_t r = 0; r < bases.size(); ++r) {
			if (bases[r] == 1) {
				++theta_0;
			}
		}
		theta_0 /= (bases.size() + 2.0);

		// initialize minimizer
		gsl_min_fminimizer* s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
		gsl_min_fminimizer_set(s, &f, logit(theta_0), -1000, 1000);

		int status;
		for (size_t iter = 0; iter < max_iter; ++iter) {
			status = gsl_min_fminimizer_iterate(s);
			if (status == GSL_FAILURE || status == GSL_EBADFUNC) break;
		}

		// check if minimizer converaged
		double theta_hat = 0.0;
		double a = gsl_min_fminimizer_x_lower(s);
		double b = gsl_min_fminimizer_x_upper(s);
		if (gsl_min_test_interval(a, b, epsabs, 0.0) == GSL_SUCCESS) {
			theta_hat = logistic(gsl_min_fminimizer_x_minimum(s));
		}

		gsl_min_fminimizer_free(s);

		return theta_hat;
	}

	/**
	 * Calculate deviance of fitted and null models.
	 *
	 * model_1: theta = theta.hat
	 * model_0: theta = 0
	 */
	double deviance_theta(double theta_hat, double phi) {
		return 2 * (lp_bases_given(theta_hat, phi) - lp_bases_given(0.0, phi));
	}

	/**
	 * Variant test with adjustment for orientation bias.
	 *
	 * Test whether variant allele frequence > 0,
	 * with consideration for orietation bias, as specified
	 * by the estimated global damage rate causing  orientation bias.
	 *
	 * Relevant reads should have been tallied using push() already.
	 *
	 * @param phi estimated global damage rate
	 * @return p-value
	 */
	double operator()(double phi) {
		double theta_hat = estimate_theta_given(phi);
		double dev = deviance_theta(theta_hat, phi);
		return 1 - gsl_cdf_chisq_P(dev, 1);
	}

};

// theta_real is in unconstrained, real-number space
// param should point to orient_bias_filter_f
inline double _nlp_bases_given(double theta_real, void* p) {
	orient_bias_filter_f* x = (orient_bias_filter_f*) p;
	return - x->lp_bases_given(logistic(theta_real), x->phi);
}

}  // namespace hts

#endif  // _HTSPAN_ORIENT_BIAS_HPP_
