#ifndef _HTSPAN_BAYES_ORIENT_BIAS_QUANT_HPP_
#define _HTSPAN_BAYES_ORIENT_BIAS_QUANT_HPP_

#include <ostream>
#include <stdexcept>
#include <cstdlib>

#include <alglib/ap.h>
#include <alglib/alglibmisc.h>
#include <alglib/specialfunctions.h>

#include <stograd/src/stograd/stograd.hpp>

#include "math.hpp"
#include "base_orient_bias_quant.hpp"
#include "freq_orient_bias_quant.hpp"

namespace hts {

using namespace std;

/**
* Convenience class to hold parameter values for Bayesian quantification.
*/
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

struct bayes_quant_model {

	// vector of consistent alt counts, by site
	vector<long int> xc_vec;

	// vector of inconsistent alt counts, by site
	vector<long int> xi_vec;

	// vector of consistent total counts, by site
	vector<long int> nc_vec;

	//vector of inconsistent total counts, by site
	vector<long int> ni_vec;

	// index of current site being examined
	size_t J;

	bayes_quant_model () : 
		xc_vec(0), xi_vec(0), nc_vec(0), ni_vec(0), J(0) {}

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

	/**
	* Total number of sites in data set.
	*/
	size_t nobs () const {
		return xc_vec.size();
	}

	/*
	* Advance to the next site for analysis.
	*/
	size_t next () {
		if (J == nobs()-1) {
			J = 0;
		} else {
			++J;
		}
		return J;
	}

	/*
	* Gradient of the objective function for optimizing
	* the parameters of the beta distribution that characterizes
	* theta (alt allele probability), with respect to alpha.
	*/
	static double dlp_xij_given_hparams_dalpha (long int xij, long int nij, double alpha_theta, double beta_theta) {
		static double alpha_cached = -1.0;
		static double beta_cached = -1.0;
		static double psi_alpha_beta = -1.0;
		static double psi_alpha = -1.0;
		if (alpha_theta != alpha_cached) {
			psi_alpha_beta = alglib::psi(alpha_theta + beta_theta);
			psi_alpha = alglib::psi(alpha_theta);
			alpha_cached = alpha_theta;
			beta_cached = beta_theta;
		} else if (beta_theta != beta_cached) {
			psi_alpha_beta = alglib::psi(alpha_theta + beta_theta);
			beta_cached = beta_theta;
		}
		return psi_alpha_beta - psi_alpha + alglib::psi(alpha_theta + xij) - alglib::psi(alpha_theta + beta_theta + nij);
	}

	/*
	* Gradient of the objective function for optimizing
	* the parameters of the beta distribution that characterizes
	* theta (alt allele probability), with respect to beta.
	*/
	static double dlp_xij_given_hparams_dbeta (long int xij, long int nij, double alpha_theta, double beta_theta) {
		static double alpha_cached = -1.0;
		static double beta_cached = -1.0;
		static double psi_alpha_beta = -1.0;
		static double psi_beta = -1.0;
		if (beta_theta != beta_cached) {
			psi_alpha_beta = alglib::psi(alpha_theta + beta_theta);
			psi_beta = alglib::psi(beta_theta);
			alpha_cached = alpha_theta;
			beta_cached = beta_theta;
		} else if (alpha_theta != alpha_cached) {
			psi_alpha_beta = alglib::psi(alpha_theta + beta_theta);
			alpha_cached = alpha_theta;
		}
		return psi_alpha_beta - psi_beta + alglib::psi(beta_theta + nij - xij) - alglib::psi(alpha_theta + beta_theta + nij);
	}

	/*
	* Evaluate the gradients of the theta objective function
	* wrt alpha and beta; x[0] == alpha and x[1] == beta
	*/
	vector<double> theta_gradient (vector<double> x) {
		vector<double> grad (2);
		grad[0] = dlp_xij_given_hparams_dalpha(xi_vec[J], ni_vec[J], x[0], x[1]);
		grad[1] = dlp_xij_given_hparams_dbeta(xi_vec[J], ni_vec[J], x[0], x[1]);
		return grad;
	}

	/*
	* The objective function used for optimizing the parameters
	* of the beta distribution that characterizes phi (damage probability).
	*/
	static double lp_xcj_given_hparams (long int xcj, long int ncj, double alpha_phi, double beta_phi,
			const double alpha_theta, const double beta_theta) {
		static double alpha_theta_cached = -1.0;
		static double beta_theta_cached = -1.0;
		static double alpha_phi_cached = -1.0;
		static double beta_phi_cached = -1.0;
		static double lbeta_ath_bth = -1.0;
		static double lbeta_aph_bph = -1.0;
		if (alpha_theta != alpha_theta_cached || beta_theta != beta_theta_cached) {
			lbeta_ath_bth = lbeta(alpha_theta, beta_theta);
			alpha_theta_cached = alpha_theta;
			beta_theta_cached = beta_theta;
		}
		if (alpha_phi != alpha_phi_cached || beta_phi != beta_phi_cached) {
			lbeta_aph_bph = lbeta(alpha_phi, beta_phi);
			alpha_phi_cached = alpha_phi;
			beta_phi_cached = beta_phi;
		}
		double *lse_array = new double[xcj+1];
		for (int k = 0; k <= xcj; ++k) {
			lse_array[k] =
				-log(ncj - k + 1) +
				lbeta(alpha_phi + xcj - k, beta_phi + ncj - xcj) -
				lbeta(xcj - k + 1,         ncj - xcj + 1) +
				lbeta(alpha_theta + k, beta_theta + ncj - k) -
				lbeta(k + 1, ncj - k + 1);
		}
		double lse = log_sum_exp(xcj+1, lse_array);
		delete[] lse_array;
		return lse - log(ncj + 1) - lbeta_ath_bth - lbeta_aph_bph;
	}

	/**
	* Functor used to numerically evaluate the gradient of 
	* the phi objective function, since the analytic gradient is
	* expensive to evaluate.
	*
	* A functor is used because stograd::finite_difference_gradient
	* expects one.
	*/
	struct phi_objective_f {
		const bayes_quant_model &m;

		const double alpha_theta;

		const double beta_theta;

		phi_objective_f(const bayes_quant_model &_m, const double a_th, const double b_th) :
			m(_m), alpha_theta(a_th), beta_theta(b_th) {}

		double operator() (const vector<double> &x) const {
			return m.lp_xcj_given_hparams(m.xc_vec[m.J], m.nc_vec[m.J],
				x[0], x[1], alpha_theta, beta_theta);
		}
	};

	/*
	* Evaluate the gradients of the phi objective function
	* wrt alpha and beta; x[0] == alpha and x[1] == beta
	*/
	vector<double> phi_gradient (vector<double> x, const double alpha_theta, const double beta_theta) {
		static phi_objective_f F (*this, alpha_theta, beta_theta);
		vector<double> grad;
		stograd::finite_difference_gradient(F, x, grad);// grad passed by reference and holds return values
		return grad;
	}
};

struct theta_hparams_optimizable {
	bayes_quant_model &m;

	vector<double> curr_params;

	theta_hparams_optimizable(bayes_quant_model &_m, vector<double> initial_params)
		: m(_m), curr_params(log_c(initial_params)) {}

	size_t nobs() const {
		return m.nobs();
	}

	size_t nparams() const {
		return 2;
	}

	double alpha() const {
		return exp(curr_params[0]);
	}

	double beta() const {
		return exp(curr_params[1]);
	}

	void accumulate (vector<double> &current_grad) {
		vector<double> new_grad = m.theta_gradient(exp_c(curr_params));
		if (isnan(new_grad[0]) || isnan(new_grad[1])) {
			throw runtime_error ("NaN returned from gradient!");
		}
		stograd::add_to (new_grad, current_grad);
		m.next();
	}

	void update (const vector<double> &delta) {
		if (isnan(delta[0]) || isnan(delta[1])) {
			throw runtime_error ("NaN returned from stograd!");
		}
		stograd::add_to(delta, curr_params);
		constrain_values(curr_params, -15, 15);
	}

};

struct phi_hparams_optimizable {
	bayes_quant_model &m;

	vector<double> curr_params;

	const double alpha_theta;
	const double beta_theta;

	phi_hparams_optimizable (bayes_quant_model &_m, vector<double> initial_params, double a_t, double b_t) :
		m(_m), curr_params(log_c(initial_params)), alpha_theta(a_t), beta_theta(b_t) {}

	size_t nobs() const {
		return m.nobs();
	}

	size_t nparams() const {
		return 2;
	}

	double alpha() const {
		return exp(curr_params[0]);
	}

	double beta() const {
		return exp(curr_params[1]);
	}

	void accumulate (vector<double> &current_grad) {
		vector<double> new_grad = m.phi_gradient(exp_c(curr_params), alpha_theta, beta_theta);
		if (isnan(new_grad[0]) || isnan(new_grad[1])) {
			throw runtime_error ("NaN returned from gradient!");
		}
		stograd::add_to(new_grad, current_grad);
		m.next();
	}

	void update (const vector<double> &delta) {
		if (isnan(delta[0]) || isnan(delta[1])) {
			throw runtime_error ("NaN returned from stograd!");
		}
		stograd::add_to(delta, curr_params);
		constrain_values(curr_params, -15, 15);
	}

};


struct bayes_orient_bias_quant_f : public base_orient_bias_quant_f {

	bayes_quant_model m;

	int n_phi_epochs;
	int n_theta_epochs;

	bayes_orient_bias_quant_f(nuc_t _ref, nuc_t _alt) :
			base_orient_bias_quant_f(_ref, _alt),
			n_phi_epochs(0), n_theta_epochs(0)
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
		m.xc_vec.push_back(xc);
		m.xi_vec.push_back(xi);
		m.nc_vec.push_back(nc);
		m.ni_vec.push_back(ni);
		read_count += success;
		return success;
	}

	/*
	* Generate a set of appropriately distributed observed variables, for testing.
	* In the Bayesian model, phi and theta are not constant but instead beta-distributed.
	* 
	* @param N Number of loci to simulate
	* @param alpha_theta Alpha parameter of the distribution of theta
	* @param beta_theta Beta parameter of the distribution of theta
	* @param alpha_phi Alpha parameter of the distribution of phi
	* @param beta_phi Beta parameter of the distribution of phi
	*/
	void simulate(size_t N, double alpha_theta, double beta_theta,
			double alpha_phi, double beta_phi, size_t rd_min, size_t rd_max, int seed = 0) {
		vector<double> theta_vec (N);
		vector<double> phi_vec (N);
		// Generate uniformly distributed random read counts
		alglib::hqrndstate rng;
		srand(seed);
		int s1 = rand();
		int s2 = rand();
		alglib::hqrndseed(s1, s2, rng);
		vector<size_t> ns_vec (N);
		for (size_t j = 0; j < N; ++j) {
			ns_vec[j] = alglib::hqrnduniformi(rng, rd_max - rd_min) + rd_min;
		}
		// Generate correctly distributed values of phi and theta
		// to pass to the simulator
		for (size_t j = 0; j < N; ++j) {
			theta_vec[j] = r_rand::rbeta(alpha_theta, beta_theta, seed);
		}
		for (size_t j = 0; j < N; ++j) {
			phi_vec[j] = r_rand::rbeta(alpha_phi, beta_phi, seed);
		}
		// Simulate as many reads as necessary
		m.reset_realloc(N);
		for (size_t j = 0; j < N; ++j) {
			simulate_orient_bias_read_counts(ns_vec[j], theta_vec[j], phi_vec[j], xc, xi, nc, ni, seed);
			m.xc_vec.push_back(xc);
			m.xi_vec.push_back(xi);
			m.nc_vec.push_back(nc);
			m.ni_vec.push_back(ni);
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
	hparams operator()(const size_t bsize, const size_t nepochs, const double learning_rate, const double eps,
			const double alpha0_theta, const double beta0_theta, const double alpha0_phi, const double beta0_phi) {
		vector<double> theta_init(2);
		theta_init[0] = alpha0_theta;
		theta_init[1] = beta0_theta;
		theta_hparams_optimizable theta_opt (m, theta_init);
		bayes_stepper_t stepper (learning_rate);
		n_theta_epochs = stograd::optimize(theta_opt, stepper, bsize, nepochs, eps);
		vector<double> phi_init(2);
		phi_init[0] = alpha0_phi;
		phi_init[1] = beta0_phi;
		phi_hparams_optimizable phi_opt (m, phi_init, theta_opt.alpha(), theta_opt.beta());
		m.J = 0;
		n_phi_epochs = stograd::optimize(phi_opt, stepper, bsize, nepochs, eps);
		hparams rtn;
		rtn.alpha_theta = theta_opt.alpha();
		rtn.beta_theta = theta_opt.beta();
		rtn.alpha_phi = phi_opt.alpha();
		rtn.beta_phi = phi_opt.beta();
		return rtn;
	}

	template <typename bayes_stepper_t>
	hparams operator() (const size_t bsize, const size_t nepochs, const double learning_rate, const double eps) {
		// obtain the frequentist estimate to use as a starting point, assuming alpha = 1
		freq_orient_bias_quant_f fobquant (r1_ref, r1_alt);
		fobquant.copy_data(m.xc_vec, m.xi_vec, m.nc_vec, m.ni_vec);
		double beta0_theta = (double) 1 / fobquant.theta_hat() - 1;
		double beta0_phi = (double) 1 / fobquant() - 1; // fobquant.operator() returns phi_hat
		return operator()<bayes_stepper_t> (bsize, nepochs, learning_rate, eps,
				1, beta0_theta, 1, beta0_phi);
	}

};// struct bayes_orient_bias_quant_f

} // namespace hts

#endif // _HTSPAN_BAYES_ORIENT_BIAS_QUANT_HPP_