#ifndef _HTSPAN_ORIENT_BIAS_QUANT_HPP_
#define _HTSPAN_ORIENT_BIAS_QUANT_HPP_

#include <ostream>
#include <vector>
#include <cassert>
#include <cfloat>
#include <utility>

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


using namespace std;

struct hparams {
	double alpha_theta;
	double beta_theta;
	double alpha_phi;
	double beta_phi;

	// Return parameter vectors for theta or phi
	vector<double> theta_pvec() {
		vector<double> rtn (2);
		rtn[0] = alpha_theta;
		rtn[1] = beta_theta;
		return rtn;
	}

	vector<double> phi_pvec() {
		vector<double> rtn (2);
		rtn[0] = alpha_phi;
		rtn[1] = beta_phi;
		return rtn;
	}
};

/**
* The base class for the frequentist and Bayesian quantification
* functors. Contains individual read push and simulation code.
*/
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

	/*
	* The total number of reads pushed (not loci).
	* size() returns the number of loci instead.
	*/
	size_t n_reads() const {
		return read_count;
	}

	/**
	 * Process a single read and update observed variables accordingly.
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

	/*
	* Simulate a correctly distributed set of observed variables for a read, for testing.
	*
	* @param N Number of reads at the locus
	* @param theta Alternate allele probability at the locus
	* @param phi Damage probability at the locus
	* @param xc returned as damage-consistent alt count at the locus
	* @param xi returned as damage-inconsistent alt count at the locus
	* @param nc returned as correctly oriented alt count at the locus
	* @param ni returned as incorrectly oriented alt count at the locus
	*/
	static void simulate_orient_bias_read_counts (const size_t N, const double theta, const double phi,
			long int &xc, long int &xi, long int &nc, long int &ni) {
		ni = N/2;
		nc = N - ni;
		xi = r_rand::rbinom(ni, theta);
		long int xr = r_rand::rbinom(nc, theta);
		long int xd = r_rand::rbinom(nc - xr, phi);
		xc = xr + xd;
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

	// count of read loci (not reads)
	size_t site_count;

	freq_orient_bias_quant_f(nuc_t _ref, nuc_t _alt) : 
			base_orient_bias_quant_f(_ref, _alt),
			site_count(0)
		{
		}

	/*
	* Returns the count of read loci (not reads).
	* n_reads() returns total reads instead.
	*/
	size_t size () const {
		return site_count;
	}

	void copy_data (vector<long int> xc_vec, vector<long int> xi_vec, vector<long int> nc_vec, vector<long int> ni_vec) {
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
	 * Process a bam pileup object containing
	 * the reads at one locus and update the observed variables accordingly.
	 *
	 * @param pile bam pileup object
	 * @param n    number of reads in pileup object
	 * @param pos  locus reference position
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
	void simulate(size_t N, double theta, double phi) {
		simulate_orient_bias_read_counts(N, theta, phi, xc, xi, nc, ni);
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

	/*
	* Process the reads at one locus, contained in a BAM pileup object.
	*
	* @param pile Already populated BAM pileup object 
	* @param n Number of reads in pileup
	* @param pos Reference position of the locus
	* @return Nubmer of successfully processed reads
	*/
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
		read_count += success;
		return success;
	}

	/*
	* Clear data vectors and reserve the requested
	* number of entries in them.
	*
	* @param N Number of entries to reserve
	*/ 
	void reset_realloc (size_t N) {
		xc_vec.clear();
		xc_vec.reserve(N);
		xi_vec.clear();
		xi_vec.reserve(N);
		nc_vec.clear();
		nc_vec.reserve(N);
		ni_vec.clear();
		ni_vec.reserve(N);
	}

	/*
	* Returns the number of processed loci (not total reads).
	* For total reads, use n_reads().
	*/
	size_t size() const {
		return xc_vec.size();
	}

	/*
	* Generate a set of appropriately distributed observed variables, for testing.
	* In the Bayesian model, phi and theta are not constant but instead beta-distributed.
	* 
	* @param ns_vec Vector of read counts at each simulated locus
	* @param alpha_theta Alpha parameter of the distribution of theta
	* @param beta_theta Beta parameter of the distribution of theta
	* @param alpha_phi Alpha parameter of the distribution of phi
	* @param beta_phi Beta parameter of the distribution of phi
	*/
	pair<double, double> simulate(const vector<size_t> &ns_vec, double alpha_theta, double beta_theta,
			double alpha_phi, double beta_phi) {
		vector<double> theta_vec (ns_vec.size());
		vector<double> phi_vec (ns_vec.size());
		double phi_sum = 0.0;
		double theta_sum = 0.0;
		// Generate correctly distributed values of phi and theta
		// to pass to the simulator
		for (size_t j = 0; j < ns_vec.size(); ++j) {
			theta_vec[j] = r_rand::rbeta(alpha_theta, beta_theta);
			theta_sum += theta_vec[j];
		}
		for (size_t j = 0; j < ns_vec.size(); ++j) {
			phi_vec[j] = r_rand::rbeta(alpha_phi, beta_phi);
			phi_sum += phi_vec[j];
		}
		// Simulate as many reads as necessary
		reset_realloc(ns_vec.size());
		for (size_t j = 0; j < ns_vec.size(); ++j) {
			simulate_orient_bias_read_counts(ns_vec[j], theta_vec[j], phi_vec[j], xc, xi, nc, ni);
			xc_vec.push_back(xc);
			xi_vec.push_back(xi);
			nc_vec.push_back(nc);
			ni_vec.push_back(ni);
		}
		return make_pair(theta_sum/ns_vec.size(), phi_sum/ns_vec.size());
	}

	void write_data (ostream &prn) {
		prn << "xc\txi\tnc\tni\n";
		for (size_t j = 0; j < size(); ++j) {
			prn << xc_vec[j] << '\t' << xi_vec[j] << '\t' << nc_vec[j] << '\t' << ni_vec[j] << '\n';
		}
	}

	/*
	* Estimate the alpha and beta parameters that characterize
	* phi, for use in Bayesian identification model.
	*
	* Uses stochastic gradient optimization method (stograd).
	*
	* @template bayes_stepper_t stograd::stepper object to use for optimization
	* @param bsize Batch size for stograd
	* @param nepochs Number of epochs for stograd
	* @param learning_rate Learning rate for stograd
	* @param eps Convergence threshold for stograd
	* @param alpha0_theta Initial estimate of alpha_theta
	* @param beta0_theta Initial estimate of beta_theta
	* @param alpha0_phi Initial estimate of alpha_phi
	* @param beta0_phi Initial estimate of beta_phi
	* @return A hts::hparams object containing estimates of the four hyperparameters
	*/
	template <typename bayes_stepper_t>
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

	/*
	* Evaluate the theta objective function (-lp_xi_given_hparams) at
	* the alpha_theta and beta_theta in the passed hparams object,
	* taking into account read sites up to maxJ [default: all].
	*/
	double eval_theta_objective_func (hparams x, int maxJ = -1) {
		theta_hparams_optimizable F (*this, x.theta_pvec());
		if (maxJ == -1) {
			F.J = size()-1;
		} else {
			F.J = maxJ;
		}
		return F(x.theta_pvec(), false);
	}

	/*
	* Evaluate the theta gradient wrt alpha_theta at
	* the alpha_theta and beta_theta in the passed hparams object,
	* taking into account read sites up to maxJ [default: all].
	*/
	double eval_theta_grad_dalpha (hparams x, int maxJ = -1) {
		theta_hparams_optimizable F(*this, x.theta_pvec());
		if (maxJ == -1) {
			F.J = size()-1;
		} else {
			F.J = maxJ;
		}
		return F.dlp_xi_given_hparams_dalpha(x.theta_pvec(), false);
	}

	/*
	* Evaluate the theta gradient wrt beta_theta at
	* the alpha_theta and beta_theta in the passed hparams object,
	* taking into account read sites up to maxJ [default: all].
	*/
	double eval_theta_grad_dbeta (hparams x, int maxJ = -1) {
		theta_hparams_optimizable F(*this, x.theta_pvec());
		if (maxJ == -1) {
			F.J = size()-1;
		} else {
			F.J = maxJ;
		}
		return F.dlp_xi_given_hparams_dbeta(x.theta_pvec(), false);
	}

	/*
	* Evaluate the phi objective function (-lp_xc_given_hparams) at
	* the point indicated by the passed hparams object,
	* taking into account read sites up to maxJ [default: all].
	*/
	double eval_phi_objective_func (hparams x, int maxJ = -1) {
		phi_hparams_optimizable F(*this, x.phi_pvec(), x.alpha_theta, x.beta_theta);
		if (maxJ == -1) {
			F.J = size()-1;
		} else {
			F.J = maxJ;
		}
		return F(x.phi_pvec(), false);
	}

	/*
	* Replace the stored data with that in the passed vectors.
	* All the arguments must be the same length.
	*/
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

	/*
	* Optimizable functor for use in stograd.
	*
	* Optimizes the alpha and beta parameters of the 
	* beta distribution that characterizes theta (alt allele probability).
	*
	* Stochastic gradient optimization takes an objective function, 
	* and its gradient evaluated either analytically or numerically
	* and returns a set of values to update each parameter being optimized.
	*
	* A subset of observations (the observed variables for each locus) are under
	* consideration at any one time; after each gradient accumulation, another
	* locus is added to the set being considered.	
	* 
	*/
	struct theta_hparams_optimizable {

		// parent object containing data
		const bayes_orient_bias_quant_f &parent;

		// index of the highest-indexed locus currently in the observed set
		size_t J;

		// vector of parameters being optimized
		// values are in unconstrained (log) space
		// [0] = alpha_theta
		// [1] = beta_theta
		vector<double> curr_params;

		/**
		* Values of alpha and beta are stored in log space,
		* so we transform them to the constrained space before
		* returning them.
		*/
		double alpha() const {
			return exp(curr_params[0]);
		}

		double beta() const {
			return exp(curr_params[1]);
		}

		theta_hparams_optimizable (bayes_orient_bias_quant_f &P, vector<double> initial_params) :
			parent(P), J(0), curr_params(log_c(initial_params)) {}

		/**
		* Total number of loci for which data is 
		* present in the parent.
		*/
		size_t nobs() const {
			return parent.xi_vec.size();
		}

		/**
		* The highest-indexed locus currently being
		* considered, indexed from one. Advanced by next().
		*/
		size_t cobs() const {
			return J+1;
		}

		/**
		* Number of parameters being optimized. Constant.
		*/
		size_t nparams() const {
			return 2;
		}

		/**
		* Add another locus (set of observed variables)
		* to the set of loci being considered.
		*/
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

		// Gradient of lp_xi_given_hparams wrt alpha
		double dlp_xi_given_hparams_dalpha (vector<double> x, bool exp_args = true) {
			static double alpha_cached = -1.0;
			static double psi_alpha_cached = 0.0;
			if (exp_args) {
				exp_elements(x);
			}
			if (x[0] != alpha_cached) {
				psi_alpha_cached = alglib::psi(x[0]);
				alpha_cached = x[0];
			}
			// xi_vec and ni_vec are class members
			double sum_dlp_xij_dalpha = 0.0;
			for (size_t j = 0; j < cobs(); ++j) {
				sum_dlp_xij_dalpha += dlp_xij_given_hparams_dalpha(parent.xi_vec[j], parent.ni_vec[j], x[0], x[1]);
			}
			return sum_dlp_xij_dalpha + cobs()*(alglib::psi(x[0] + x[1]) - psi_alpha_cached);
		}

		// Gradient of lp_xi_given_hparams wrt beta
		double dlp_xi_given_hparams_dbeta (vector<double> x, bool exp_args = true) {
			static double beta_cached = -1.0;
			static double psi_beta_cached = 0.0;
			if (exp_args) {
				exp_elements(x);
			}
			if (x[1] != beta_cached) {
				psi_beta_cached = alglib::psi(x[1]);
				beta_cached = x[1];
			}
			double sum_dlp_xij_dbeta = 0.0;
			for (size_t j = 0; j < cobs(); ++j) {
				sum_dlp_xij_dbeta += dlp_xij_given_hparams_dbeta(parent.xi_vec[j], parent.ni_vec[j], x[0], x[1]);
			}
			return sum_dlp_xij_dbeta + cobs()*(alglib::psi(x[0] + x[1]) - psi_beta_cached);
		}

		double lp_xi_given_hparams (vector<double> x) const {
			// xi_vec and ni_vec are class members
			// x[0] is alpha_theta, x[1] is beta_theta
			double sum_lp_xij = 0.0;
			for (size_t j = 0; j < cobs(); ++j) {
				sum_lp_xij += lp_xij_given_hparams(parent.xi_vec[j], parent.ni_vec[j], x[0], x[1]);
			}
			return sum_lp_xij - lbeta(x[0], x[1])*cobs();
		}

		/**
		* Objective function being minimized.
		*
		* Maximizing lp_xi_given_hparams, so we minimize
		* -lp_xi_given_hparams.
		*/
		double operator() (vector<double> x, bool exp_args = true) const {
			if (exp_args) {
				exp_elements(x);
			}
			return -lp_xi_given_hparams(x);
		}

		/**
		* Evaluate the gradient for the current observed set of loci
		* and add it to the accumulated gradient tracked by the optimizer,
		* then add the next locus to the observed set.
		*/
		void accumulate (vector<double> &current_grad) {
			vector<double> new_grad (2);
			new_grad[0] = dlp_xi_given_hparams_dalpha(curr_params);
			new_grad[1] = dlp_xi_given_hparams_dbeta(curr_params);
			stograd::add_to(new_grad, current_grad);
			next();
		}

		/**
		* Update the parameter values with the delta returned 
		* by the optimizer.
		*/
		void update (vector<double> &delta) {
			stograd::subtract_from(delta, curr_params);
		}
	}; // struct theta_hparams_optimizable

	/*
	* Optimizable functor for use in stograd.
	*
	* Optimizes the alpha and beta parameters of
	* the beta distribution that characterizes phi (damage probability).
	*
	* See theta_hparams_optimizable description above for details on
	* the stograd optimization method.
	*/
	struct phi_hparams_optimizable {

		// parent object containing data
		const bayes_orient_bias_quant_f &parent;

		// fixed parameters for objective function (previously estimated)
		// values are in real space
		const double alpha_theta;
		const double beta_theta;

		// vector of parameters being optimized
		// values are in unconstrained (log) space 
		// [0] = alpha_phi
		// [1] = beta_phi
		vector<double> curr_params;

		/**
		* Values of alpha and beta are stored in log space,
		* so we transform them to the constrained space before
		* returning them.
		*/
		double alpha() const {
			return exp(curr_params[0]);
		}

		double beta() const {
			return exp(curr_params[1]);
		}

		// index of record (in parent object) currently being analyzed
		size_t J;

		phi_hparams_optimizable(const bayes_orient_bias_quant_f &P, vector<double> initial_params,
				const double a_theta, const double b_theta) :
				parent(P), alpha_theta(a_theta), beta_theta(b_theta), curr_params(log_c(initial_params)), J(0) {}

		/**
		* Total number of loci for which data is present in the parent.
		*/
		size_t nobs() const {
			return parent.size();
		}

		/*
		* The highest-indexed locus currently being
		* considered. Advanced by next().
		*/
		size_t cobs() const {
			return J+1;
		}

		/*
		* Number of parameters being optimized. Constant.
		*/
		size_t nparams() const {
			return 2;
		}

		/*
		* Add another locus (set of observed variables)
		* to the set of loci being considered.
		*/
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

		double lp_xc_given_hparams (vector<double> x) const {
			double sum_lp_xcj = 0.0;
			double sum_log_ncp = 0.0;
			for (size_t j = 0; j < cobs(); ++j) {
				sum_lp_xcj += lp_xcj_given_hparams(parent.xc_vec[j], parent.nc_vec[j], alpha_theta, beta_theta, x[0], x[1]);
				sum_log_ncp += log(parent.nc_vec[j] + 1);
			}
			return sum_lp_xcj - sum_log_ncp - cobs()*(lbeta(alpha_theta, beta_theta) + lbeta(x[0], x[1]));
		}

		/**
		* Objective function being optimized.
		*
		* We want to maximize lp_xc_given_hparams, so
		* we minimize its negative.
		*/
		double operator() (vector<double> x, bool exp_args = true) const {
			if (exp_args) {
				exp_elements(x);
			}
			return -lp_xc_given_hparams(x);
		}

		/**
		* Evaluate the gradient for the current observed set of loci
		* and add it to the accumulated gradient tracked by the optimizer,
		* then add the next locus to the observed set.
		*/
		void accumulate (vector<double> &current_grad) {
			vector<double> new_grad;
			stograd::finite_difference_gradient(*this, curr_params, new_grad);
			stograd::add_to(new_grad, current_grad);
			next();
		}

		/**
		* Update the parameter values with the delta returned 
		* by the optimizer.
		*/
		void update (vector<double> &delta) {
			stograd::subtract_from(delta, curr_params);
		}
	};// struct phi_hparams_optimizable

};// struct bayes_orient_bias_quant_f

}  // namespace hts

#endif  // _HTSPAN_ORIENT_BIAS_QUANT_HPP_
