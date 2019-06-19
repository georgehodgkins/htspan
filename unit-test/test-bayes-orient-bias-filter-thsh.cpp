#define BOOST_TEST_MODULE bayes_orient_bias_filter_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#define COMMON_ALPHA 0.1
#define COMMON_BETA 0.1
#define COMMON_PRIOR 0.5
#define COMMON_PHI 0.01
#define COMMON_THETA 0.1

#include <cmath>

#include <vector>
#include <iostream>

#include "htspan/bayes_orient_bias_filter.hpp"
#include "htspan/de_integrator.hpp"
#define INTEGRATOR math::tanh_sinh<math::numeric_functor>
#include "htspan/nucleotide.hpp"

#include "test.hpp"

using namespace std;

void common_model_test(const char TSVNAME[], const double ALPHA_PHI, const double BETA_PHI,
		const double PRIOR_ALT, const double PHI_T, const double THETA_T, 
		const double PHI_INT_STD, const double THETA_INT_STD, const double EV_NULL_STD,
		const double EV_ALT_STD, const double LPOSTERIOR_STD) {
	BOOST_TEST_MESSAGE("Running model test for dataset " << TSVNAME << ":");
	BOOST_TEST_CHECKPOINT("Initializing filter and reading in data");
	hts::orient_bias_data data(nuc_G, nuc_T, 200);
	data.read(TSVNAME);
	BOOST_TEST_CHECKPOINT("Evaluating phi integrand");
	hts::bayes_orient_bias_filter_f bobfilter (data);
	hts::bayes_orient_bias_filter_f::lp_bases_theta_phi_f p_f (bobfilter, ALPHA_PHI, BETA_PHI, THETA_T);
	double phi_int = p_f(PHI_T);
	BOOST_CHECK_MESSAGE(test_val(phi_int, PHI_INT_STD),
		"Phi integrand: got: " << phi_int << ", expected: " << PHI_INT_STD);
	BOOST_TEST_CHECKPOINT("Evaluating theta integrand");
	hts::bayes_orient_bias_filter_f::lp_bases_theta_f<INTEGRATOR> t_f (bobfilter, ALPHA_PHI, BETA_PHI);
	double theta_int = t_f(THETA_T);
	BOOST_CHECK_MESSAGE(test_val(theta_int, THETA_INT_STD),
		"Theta integrand: got: " << theta_int << ", expected: " << THETA_INT_STD);
	BOOST_TEST_CHECKPOINT("Evaluating evidence");
	hts::evidences ev = bobfilter.model_evidence<INTEGRATOR>(ALPHA_PHI, BETA_PHI);
	BOOST_CHECK_MESSAGE(test_val(log(ev.null), EV_NULL_STD),
		"Evidence for null model: got: " << log(ev.null) << ", expected: " << EV_NULL_STD);
	BOOST_CHECK_MESSAGE(test_val(log(ev.alt), EV_ALT_STD),
		"Evidence for alternate model: got: " << log(ev.alt) << ", expected: " << EV_ALT_STD);
	BOOST_TEST_CHECKPOINT("Evaluating posterior probability");
	double lposterior = bobfilter.operator()<INTEGRATOR>(PRIOR_ALT, ALPHA_PHI, BETA_PHI);
	BOOST_CHECK_MESSAGE(test_val(lposterior, LPOSTERIOR_STD),
		"Log posterior probability: got: " << lposterior << ", expected: " << LPOSTERIOR_STD);
} 

BOOST_AUTO_TEST_SUITE(test)

BOOST_AUTO_TEST_CASE(hs_hd) {
	const char TSVNAME[] = "../sim-data/obrs_high-signal_high-damage_data.tsv";
	const double ALPHA_PHI = COMMON_ALPHA;
	const double BETA_PHI = COMMON_BETA;
	const double PRIOR_ALT = COMMON_PRIOR;
	const double PHI_T = COMMON_PHI;
	const double THETA_T = COMMON_THETA;
	const double PHI_INT_STD = 3.9844623e-32;
	const double THETA_INT_STD = 3.4907991e-33;
	const double EV_NULL_STD = -112.5565153;
	const double EV_ALT_STD = -77.3509179;
	const double LPOSTERIOR_STD = 0.0;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT, PHI_T, THETA_T,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_CASE(hs_ld) {
	const char TSVNAME[] = "../sim-data/obrs_high-signal_low-damage_data.tsv";
	const double ALPHA_PHI = COMMON_ALPHA;
	const double BETA_PHI = COMMON_BETA;
	const double PRIOR_ALT = COMMON_PRIOR;
	const double PHI_T = COMMON_PHI;
	const double THETA_T = COMMON_THETA;
	const double PHI_INT_STD = 1.8932271e-30;
	const double THETA_INT_STD = 5.2998802e-32;
	const double EV_NULL_STD = -131.9902998;
	const double EV_ALT_STD = -74.8719863;
	const double LPOSTERIOR_STD = 0.0;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT, PHI_T, THETA_T,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_CASE(ls_hd) {
	const char TSVNAME[] = "../sim-data/obrs_low-signal_high-damage_data.tsv";
	const double ALPHA_PHI = COMMON_ALPHA;
	const double BETA_PHI = COMMON_BETA;
	const double PRIOR_ALT = COMMON_PRIOR;
	const double PHI_T = COMMON_PHI;
	const double THETA_T = COMMON_THETA;
	const double PHI_INT_STD = 9.0136431e-22;
	const double THETA_INT_STD = 3.4159817e-23;
	const double EV_NULL_STD = -43.131369;
	const double EV_ALT_STD = -47.6025126;
	const double LPOSTERIOR_STD = -4.4837466;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT, PHI_T, THETA_T,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_CASE(ms_hd) {
	const char TSVNAME[] = "../sim-data/obrs_med-signal_high-damage_data.tsv";
	const double ALPHA_PHI = COMMON_ALPHA;
	const double BETA_PHI = COMMON_BETA;
	const double PRIOR_ALT = COMMON_PRIOR;
	const double PHI_T = COMMON_PHI;
	const double THETA_T = COMMON_THETA;
	const double PHI_INT_STD = 3.669267e-30;
	const double THETA_INT_STD = 5.0788968e-31;
	const double EV_NULL_STD = -83.1150929;
	const double EV_ALT_STD = -71.0309627;
	const double LPOSTERIOR_STD = -5.6484295e-06;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT, PHI_T, THETA_T,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_CASE(ls_ld) {
	const char TSVNAME[] = "../sim-data/obrs_low-signal_low-damage_data.tsv";
	const double ALPHA_PHI = COMMON_ALPHA;
	const double BETA_PHI = COMMON_BETA;
	const double PRIOR_ALT = COMMON_PRIOR;
	const double PHI_T = COMMON_PHI;
	const double THETA_T = COMMON_THETA;
	const double PHI_INT_STD = 8.3875786e-11;
	const double THETA_INT_STD = 5.4264592e-12;
	const double EV_NULL_STD = -8.3270583;
	const double EV_ALT_STD = -12.7868176;
	const double LPOSTERIOR_STD = -4.471258;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT, PHI_T, THETA_T,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_CASE(ms_ld) {
	const char TSVNAME[] = "../sim-data/obrs_med-signal_low-damage_data.tsv";
	const double ALPHA_PHI = COMMON_ALPHA;
	const double BETA_PHI = COMMON_BETA;
	const double PRIOR_ALT = COMMON_PRIOR;
	const double PHI_T = COMMON_PHI;
	const double THETA_T = COMMON_THETA;
	const double PHI_INT_STD = 2.3856537e-17;
	const double THETA_INT_STD = 7.7776958e-19;
	const double EV_NULL_STD = -53.4426378;
	const double EV_ALT_STD = -39.4986378;
	const double LPOSTERIOR_STD = -8.7942247e-07;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT, PHI_T, THETA_T,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_CASE(ns_hd) {
	const char TSVNAME[] = "../sim-data/obrs_no-signal_high-damage_data.tsv";
	const double ALPHA_PHI = COMMON_ALPHA;
	const double BETA_PHI = COMMON_BETA;
	const double PRIOR_ALT = COMMON_PRIOR;
	const double PHI_T = COMMON_PHI;
	const double THETA_T = COMMON_THETA;
	const double PHI_INT_STD = 2.4968668e-15;
	const double THETA_INT_STD = 8.2723427e-17;
	const double EV_NULL_STD = -25.8363411;
	const double EV_ALT_STD = -30.264146;
	const double LPOSTERIOR_STD = -4.4396749;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT, PHI_T, THETA_T,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_CASE(ns_ld) {
	const char TSVNAME[] = "../sim-data/obrs_no-signal_low-damage_data.tsv";
	const double ALPHA_PHI = COMMON_ALPHA;
	const double BETA_PHI = COMMON_BETA;
	const double PRIOR_ALT = COMMON_PRIOR;
	const double PHI_T = COMMON_PHI;
	const double THETA_T = COMMON_THETA;
	const double PHI_INT_STD = 8.9650484e-11;
	const double THETA_INT_STD = 5.3068831e-12;
	const double EV_NULL_STD = -8.2537219;
	const double EV_ALT_STD = -12.7709305;
	const double LPOSTERIOR_STD = -4.5280689;
	common_model_test(TSVNAME, ALPHA_PHI, BETA_PHI, PRIOR_ALT, PHI_T, THETA_T,
			PHI_INT_STD, THETA_INT_STD, EV_NULL_STD, EV_ALT_STD, LPOSTERIOR_STD);
}

BOOST_AUTO_TEST_SUITE_END()