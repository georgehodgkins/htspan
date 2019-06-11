#ifndef _HTSPAN_ORIENT_BIAS_HPP_
#define _HTSPAN_ORIENT_BIAS_HPP_

#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <stdexcept>

#include <stdint.h>
#include <string.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_cdf.h>

#include <gsl/gsl_randist.h>//for simulation
#include "io/orient_bias_stats.hpp"

#include "fetcher.hpp"
#include "math.hpp"

namespace hts {

using namespace std;

inline double _nlp_bases_given_theta(double theta_real, void* p);
inline double _nlp_bases_given_phi(double phi_real, void* p);

struct theta_and_phi {
	double theta;
	double phi;
};

struct orient_bias_filter_f {

	/// reference base (0), alternative base (1), or other base (2)
	vector<int8_t> bases;
	/// base error probability
	vector<double> errors;
	/// whether variant has an orientation consistent with artifact damage
	/// (false if base is not variant)
	vector<bool> orients;

	// global damage probability
	double phi_t;
	
	// alternate allele frequency (temporary)
	double theta_t;
	
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


	//upper and lower bounds for estimation
	const double minimizer_lb;//default -15.0
	const double minimizer_ub;//default 15.0
	//maximum error threshold for estimation
	const double epsabs;//default .001
	//maximum number of iterations for estimation
	const size_t max_minimizer_iter;//default 100

	orient_bias_filter_f(nuc_t _ref, nuc_t _alt, size_t n, 
		double lb = -15.0, double ub = 15.0, double eps = 1e-6, size_t max_iter = 100)
	: phi_t(0.0),
		theta_t(0.0),	
		r1_ref(_ref),
		r1_alt(_alt),
		r2_ref(nuc_complement(_ref)),
		r2_alt(nuc_complement(_alt)),
		minimizer_lb(lb),
		minimizer_ub(ub),
		epsabs(eps),
		max_minimizer_iter(max_iter)
	{
		reserve(n);
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
			members.push_back(2);
			if (onuc == r2_ref || onuc == r2_alt) {
				// e.g. C or A on read 2
				orients.push_back(true);
			} else {
				orients.push_back(false);
			}
		} else {
			// flag is malformed
			members.push_back(0);
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
	* Estimate initial value of theta from reads.
	*/
	double estimate_initial_theta() {
		double theta_0 = 1.0;
		for (size_t r = 0; r < bases.size(); ++r) {
			if (bases[r] == 1) {
				++theta_0;
			}
		}
		theta_0 /= (bases.size() + 2.0);
		return theta_0;
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
	double estimate_theta_given(double phi, double theta_0) {
		// define objective function to minimize
		phi_t = phi;
		gsl_function f;
		f.function = &_nlp_bases_given_theta;
		f.params = this;
		//call common minimization code
		return minimize_log_function(&f, theta_0);
	}
	
	/**
	*  Find value of phi that maximizes the log probability of observed bases,
	* given a value of theta.
	*/ 
	double estimate_phi_given(double theta, double phi_0) {
		// define objective function to minimize
		theta_t = theta;
		gsl_function f;
		f.function = &_nlp_bases_given_phi;
		f.params = this;
		//call common minimization code
		return minimize_log_function(&f, phi_0);
	}

	/**
	* Minimization code common to the phi and theta estimation functions, using Brent method.
	*
	* If the function value at x_0 is not less than the values at the endpoints,
	* the function will return the smaller of the endpoints and will not perform a minimization.
	*
	* epsabs, minimizer_ub, minimizer_lb, and max_minimizer_iter
	* are const class members set in constructor
	* 
	* @param f_ptr Pointer to a gsl_function struct initialized by the caller
	* @param x_0 An initial guess for the minimizer
	* @return Minimum from the GSL minimizer, if possible (see above)
	*/

	double minimize_log_function (gsl_function *f_ptr, double x_0) {
		double lx_0 = logit(x_0); // transform the guess into unconstrained space
		// if the function has already been minimized to an endpoint or
		// starts outside the bounds simply return the guess
		// generally lx_0 will only be OOB if it's -inf (x_0 == 0)
		if (lx_0 <= minimizer_lb || lx_0 >= minimizer_ub) {
			return x_0;
		}
		// if the function at the initial guess is not lower than
		// the value at both endpoints, minimizer will not work
		double f_x0 = f_ptr->function(lx_0, this);
		double f_lb = f_ptr->function(minimizer_lb, this);
		double f_ub = f_ptr->function(minimizer_ub, this);
		if (f_x0 > f_lb || f_x0 > f_ub) {
			return (f_lb < f_ub) ? logistic(minimizer_lb) : logistic(minimizer_ub);
		}
		// initialize minimizer
		gsl_min_fminimizer* s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
		gsl_min_fminimizer_set(s, f_ptr, lx_0, minimizer_lb, minimizer_ub);
		// iterate on the minimizer
		int status = GSL_CONTINUE;
		for (size_t iter = 0; iter < max_minimizer_iter; ++iter) {
			status = gsl_min_fminimizer_iterate(s);
			if (status == GSL_FAILURE || status == GSL_EBADFUNC) break;
			// test for success
			double a = gsl_min_fminimizer_x_lower(s);
			double b = gsl_min_fminimizer_x_upper(s);
			status = gsl_min_test_interval(a, b, epsabs, 0.0);
			if (status == GSL_SUCCESS) break;
		}
		// store minimizer result and then deallocate
		double xmin = logistic(gsl_min_fminimizer_x_minimum(s));
		gsl_min_fminimizer_free(s);
		switch (status) {
			case GSL_SUCCESS:
				return xmin;
			case GSL_FAILURE:
			case GSL_EBADFUNC:
				throw runtime_error("Unknown error in minimizer");
			default:
				throw runtime_error("Minimizer did not converge");
		}
		//this point should never be reached but this makes the compiler happy
		return GSL_NAN;
	}

	/**
	* Estimate theta and phi jointly using coordinate ascent.
	*
	* @param theta initial value of theta
	* @param phi initial value of phi
	* @param eps numerical tolerance
	* @param max_iter maximum number of iterations
	* @return a struct containing the estimated values of theta and phi
	*/
	theta_and_phi estimate_theta_phi(double theta, double phi, double eps = 1e-6, int max_iter = 100) {
		double theta_old = theta;
		double phi_old = phi;
		for (int n = 1; n <= max_iter; ++n) {
			theta = estimate_theta_given(phi, theta);
			phi = estimate_phi_given(theta, phi);
			if (abs(theta_old - theta) < eps && abs(phi_old - phi) < eps) {
				break;
			}
			theta_old = theta;
			phi_old = phi;
		}
		theta_and_phi rtn;
		rtn.theta = theta;
		rtn.phi = phi;
		return rtn;
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
	 * @param known_phi consider phi known/fixed (true) or variable (false)
	 * @return p-value
	 */
	double operator()(double phi, bool known_phi=true) {
		volatile double theta_0 = estimate_initial_theta();
		volatile double theta_hat;
		if (known_phi) {
			theta_hat = estimate_theta_given(phi, theta_0);
		} else {
			theta_and_phi theta_phi = estimate_theta_phi(theta_0, phi);
			theta_hat = theta_phi.theta;
			phi = theta_phi.phi;
		}
		volatile double dev = deviance_theta(theta_hat, phi);
		return 1 - gsl_cdf_chisq_P(dev, 1);
	}
	
	/**
	* C++ version of the read simulation code in R,
	* populating class members without real BAM data.
	* Used to make unit tests consistent with R model.
	*
	* Note that this will clear any reads previously pushed.
	*
	* @param theta Target theta value
	* @param phi Target phi value
	* @param err_mean Mean value for normal distribution of error values
	* @param err_sd Std dev for normal distribution of error values
	* @return none (populates class members with simulated values)
	*/
	void simulate(double theta, double phi, double err_mean, double err_sd) {
		size_t count = bases.capacity();//set in constructor
		clear();
		reserve(count);//in case clearing resets allocation, which is implementation dependent
		//NB: R docs are unclear about what RNG is used for random sampling, so this may be a point of divergence
		gsl_rng *rng = gsl_rng_alloc (gsl_rng_mt19937);
		//generate normally distributed error values with the input parameters
		for (size_t n = 0; n < count; ++n) {
			double z = gsl_ran_gaussian(rng, err_sd) + err_mean;
			if (z < 0) z = 0;
			errors[n] = anti_phred(z);
		}
		//generate bases where prob of alt is theta and prob of ref is 1-theta
		for (size_t n = 0; n < count; ++n) {
			double z = gsl_ran_flat(rng, 0, 1);
			bases[n] = (z < theta) ? 0 : 1;//0=ref, 1=alt
			//generate orients with equal probability
			z = gsl_ran_flat(rng, 0, 1);
			orients[n] = (z < 0.5);
			//for oriented reference bases, prob of damage (ref>alt) is phi
			if (bases[n] == 0 && orients[n]) {
				z = gsl_ran_flat(rng, 0, 1);
				if (z < phi) {
					bases[n] = 1;
				}
			}
		}
	}

	/**
	* Reserve space in vector members for a given number of reads.
	*
	* @param n Number of fields to reserve in each member
	* @param info_fields Whether to reserve space in informational vectors as well [false]
	*/
	void reserve (size_t n, bool info_fields = false) {
		orients.reserve(n);
		bases.reserve(n);
		errors.reserve(n);
		if (info_fields) {
			cnucs.reserve(n);
			quals.reserve(n);
			strands.reserve(n);
			members.reserve(n);
		}
	}

	/**
	* Clear member vectors.
	*/
	void clear () {
		orients.clear();
		bases.clear();
		errors.clear();
		cnucs.clear();
		quals.clear();
		strands.clear();
		members.clear();
	}

	/**
	* Reads previously generated simulation data from a .tsv file,
	* and populates class members with that data.
	*
	* Note that this will clear any reads previously pushed.
	*
	* @param path to the .tsv file
	* @return none (populates orients, bases, and errors)
	*/
	void read (const char* path) {
		clear();
		orient_bias_stats::reader stat_reader (path);
		orient_bias_stats::record rec;
		while (stat_reader.next(rec)) {
			bases.push_back(rec.base);
			errors.push_back(rec.error);
			orients.push_back(rec.orient);
		}
		stat_reader.close();
	}
};


// theta_real is in unconstrained, real-number space
inline double _nlp_bases_given_theta(double theta_real, void* p) {
	orient_bias_filter_f* x = (orient_bias_filter_f*) p;
	return - x->lp_bases_given(logistic(theta_real), x->phi_t);
}


// phi_real is in unconstrained, real-number space
inline double _nlp_bases_given_phi(double phi_real, void* p) {
	orient_bias_filter_f* x = (orient_bias_filter_f*) p;
	return - x->lp_bases_given(x->theta_t, logistic(phi_real));
}

}  // namespace hts

#endif  // _HTSPAN_ORIENT_BIAS_HPP_
