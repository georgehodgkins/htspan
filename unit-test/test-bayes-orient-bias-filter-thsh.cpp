#define BOOST_TEST_MODULE bayes_orient_bias_filter_module
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#define TEST_EPS 1e-3
#define INTEGRATOR math::tanh_sinh<math::numeric_functor>

#include <cmath>

#include <vector>
#include <iostream>

#include "htspan/bayes_orient_bias_filter.hpp"
#include "htspan/de_integrator.hpp"
#include "htspan/nucleotide.hpp"

using namespace std;

void common_model_test(const char TSVNAME[], const double ALPHA_PHI, const double BETA_PHI,
		const double PRIOR_ALT, const double PHI_INT_STD, const double THETA_INT_STD,
		const double EV_NULL_STD, const double EV_ALT_STD, const double LPOSTERIOR_STD) {
	BOOST_TEST_MESSAGE("Running model test for dataset " << TSVNAME << ":");
	BOOST_TEST_CHECKPOINT("Initializing filter and reading in data");
	hts::orient_bias_data data(nuc_G, nuc_T, 200);
	data.read(TSVNAME);
	BOOST_TEST_CHECKPOINT("Evaluating phi integrand");
	hts::bayes_orient_bias_filter_f bobfilter (data);
	// TODO: evaluate integrand fcns at fixed points, and update terminology
	double phi_hat = bobfilter.estimate_phi_given(0, 0.5); // fixed params same in model_evidence method
	double theta_hat = bobfilter.estimate_theta_given(phi_hat, 0.5); // ^^^
	hts::bayes_orient_bias_filter_f::lp_bases_theta_phi_f p_f (bobfilter, ALPHA_PHI, BETA_PHI, theta_hat);
	double phi_int = p_f(phi_hat);
	BOOST_CHECK_MESSAGE(abs(phi_int - PHI_INT_STD) < TEST_EPS,
		"Phi integrand: got: " << phi_int << ", expected: " << PHI_INT_STD);
	BOOST_TEST_CHECKPOINT("Evaluating theta integrand");
	hts::bayes_orient_bias_filter_f::lp_bases_theta_f<INTEGRATOR> t_f (bobfilter, ALPHA_PHI, BETA_PHI);
	double theta_int = t_f(theta_hat);
	BOOST_CHECK_MESSAGE(abs(theta_int - THETA_INT_STD) < TEST_EPS,
		"Theta integrand: got: " << theta_int << ", expected: " << THETA_INT_STD);
	BOOST_TEST_CHECKPOINT("Evaluating evidence");
	hts::evidences ev = bobfilter.model_evidence<INTEGRATOR>(ALPHA_PHI, BETA_PHI);
	BOOST_CHECK_MESSAGE(abs(ev.null - EV_NULL_STD) < TEST_EPS,
		"Evidence for null model: got: " << ev.null << ", expected: " << EV_NULL_STD);
	BOOST_CHECK_MESSAGE(abs(ev.alt - EV_ALT_STD) < TEST_EPS,
		"Evidence for alternate model: got: " << ev.alt << ", expected: " << EV_ALT_STD);
	BOOST_TEST_CHECKPOINT("Evaluating posterior probability");
	double lposterior = bobfilter.operator()<INTEGRATOR>(PRIOR_ALT, ALPHA_PHI, BETA_PHI);
	BOOST_CHECK_MESSAGE(abs(lposterior - LPOSTERIOR_STD) < TEST_EPS,
		"Log posterior probability: got: " << lposterior << ", expected: " << LPOSTERIOR_STD);
} 

BOOST_AUTO_TEST_SUITE(bayes_orient_bias_filter_thsh)

BOOST_AUTO_TEST_CASE(hs_hd) {
	const char TSVNAME[] = "../sim-data/obrs_high-signal_high-damage_data.tsv";
	const double ALPHA_PHI = .1;
	const double BETA_PHI = .1;
	const double PRIOR_ALT = .5;
	const double PHI_INT_STD = 9.0183999e-33;
	const double THETA_INT_STD = 2.8040116e-33;
	const double EV_NULL_STD = 1.5655737e-49;
	const double EV_ALT_STD = 3.0031175e-34;
	const double LPOSTERIOR_STD = 0.0;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_CASE(hs_ld) {
	const char TSVNAME[] = "../sim-data/obrs_high-signal_low-damage_data.tsv";
	const double ALPHA_PHI = .1;
	const double BETA_PHI = .1;
	const double PRIOR_ALT = .5;
	const double PHI_INT_STD = 2.4286403e-32;
	const double THETA_INT_STD = 4.693075e-32;
	const double EV_NULL_STD = 4.7640676e-58;
	const double EV_ALT_STD = 3.4185215e-33;
	const double LPOSTERIOR_STD = 0.0;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_CASE(ls_hd) {
	const char TSVNAME[] = "../sim-data/obrs_low-signal_high-damage_data.tsv";
	const double ALPHA_PHI = .1;
	const double BETA_PHI = .1;
	const double PRIOR_ALT = .5;
	const double PHI_INT_STD = 2.6444629e-18;
	const double THETA_INT_STD = 2.1192642e-19;
	const double EV_NULL_STD = 2.1193205e-19;
	const double EV_ALT_STD = 2.4202584e-21;
	const double LPOSTERIOR_STD = -4.4837466;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_CASE(ms_hd) {
	const char TSVNAME[] = "../sim-data/obrs_med-signal_high-damage_data.tsv";
	const double ALPHA_PHI = .1;
	const double BETA_PHI = .1;
	const double PRIOR_ALT = .5;
	const double PHI_INT_STD = 1.9561488e-29;
	const double THETA_INT_STD = 2.7535448e-30;
	const double EV_NULL_STD = 9.7007392e-37;
	const double EV_ALT_STD = 1.7065937e-31;
	const double LPOSTERIOR_STD = -5.6842537e-06;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_CASE(ls_ld) {
	const char TSVNAME[] = "../sim-data/obrs_low-signal_low-damage_data.tsv";
	const double ALPHA_PHI = .1;
	const double BETA_PHI = .1;
	const double PRIOR_ALT = .5;
	const double PHI_INT_STD = 0.0100674;
	const double THETA_INT_STD = 0.0002444;
	const double EV_NULL_STD = 0.0002444;
	const double EV_ALT_STD = 2.8219751e-06;
	const double LPOSTERIOR_STD = -4.4729571;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_CASE(ms_ld) {
	const char TSVNAME[] = "../sim-data/obrs_med-signal_low-damage_data.tsv";
	const double ALPHA_PHI = .1;
	const double BETA_PHI = .1;
	const double PRIOR_ALT = .5;
	const double PHI_INT_STD = 8.6220378e-16;
	const double THETA_INT_STD = 1.8321337e-16;
	const double EV_NULL_STD = 6.1727405e-24;
	const double EV_ALT_STD = 7.3143099e-18;
	const double LPOSTERIOR_STD = -8.4392622e-07;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_CASE(ns_hd) {
	const char TSVNAME[] = "../sim-data/obrs_no-signal_high-damage_data.tsv";
	const double ALPHA_PHI = .1;
	const double BETA_PHI = .1;
	const double PRIOR_ALT = .5;
	const double PHI_INT_STD = 1.029417e-10;
	const double THETA_INT_STD = 6.4375412e-12;
	const double EV_NULL_STD = 6.4377133e-12;
	const double EV_ALT_STD = 7.6568799e-14;
	const double LPOSTERIOR_STD = -4.4435626;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_CASE(ns_ld) {
	const char TSVNAME[] = "../sim-data/obrs_no-signal_low-damage_data.tsv";
	const double ALPHA_PHI = .1;
	const double BETA_PHI = .1;
	const double PRIOR_ALT = .5;
	const double PHI_INT_STD = 0.0102652;
	const double THETA_INT_STD = 0.0002632;
	const double EV_NULL_STD = 0.0002632;
	const double EV_ALT_STD = 2.869548e-06;
	const double LPOSTERIOR_STD = -4.5294511;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_SUITE_END()