#define BOOST_TEST_MODULE orient_bias_quant_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstddef>

#include <alglib/alglibmisc.h>

#include "../htspan/nucleotide.hpp"
#include "../htspan/orient_bias_quant.hpp"

#include "test.hpp"

#define FREQ_N_READS 10000

#define BAYES_LEARN_RATE 1e-2
#define BAYES_EPS 1e-3
#define BAYES_BSIZE 100
#define BAYES_NEPOCHS 2
#define BAYES_INIT_ALPHA .007
#define BAYES_INIT_BETA .12
#define BAYES_N_LOCS 1000
#define BAYES_MIN_RCOUNT 10
#define BAYES_MAX_RCOUNT 20

using namespace std;

void freq_obquant_sim_test (const size_t N, const double THETA, const double PHI) {
	hts::freq_orient_bias_quant_f fobquant (nuc_T, nuc_C);
	fobquant.simulate(N, THETA, PHI);
	double phi_test = fobquant();
	BOOST_CHECK_MESSAGE(test_val(phi_test, PHI),
		"Phi value returned by estimator was not within tolerance. Got: " << phi_test << " Exp: " << PHI);
}

void bayes_obquant_sim_test (const size_t N, const size_t LOC_RD_MIN, const size_t LOC_RD_MAX,
		const double ALPHA_THETA, const double BETA_THETA, const double ALPHA_PHI, const double BETA_PHI) {
	// generate random read counts within the requested range
	alglib::hqrndstate rng;
	alglib::hqrndrandomize(rng);
	vector<size_t> ns_vec (N);
	for (size_t j = 0; j < N; ++j) {
		ns_vec[j] = alglib::hqrnduniformi(rng, LOC_RD_MAX - LOC_RD_MIN) + LOC_RD_MIN;
	}
	// do the test
	hts::hparams expected;
	expected.alpha_phi = ALPHA_PHI;
	expected.beta_phi = BETA_PHI;
	expected.alpha_theta = ALPHA_THETA;
	expected.beta_theta = BETA_THETA;
	hts::bayes_orient_bias_quant_f bobquant (nuc_T, nuc_C);
	bobquant.simulate(ns_vec, ALPHA_THETA, BETA_THETA, ALPHA_PHI, BETA_PHI);
	hts::hparams alb_test = bobquant(BAYES_BSIZE, BAYES_NEPOCHS, BAYES_LEARN_RATE, BAYES_EPS,
		BAYES_INIT_ALPHA, BAYES_INIT_BETA, BAYES_INIT_ALPHA, BAYES_INIT_BETA);
	BOOST_CHECK_MESSAGE(test_val(alb_test.alpha_phi, ALPHA_PHI), 
		"Value of alpha_phi was not within tolerance. Got: " << alb_test.alpha_phi << " Exp: " << ALPHA_PHI);
	BOOST_CHECK_MESSAGE(test_val(alb_test.beta_phi, BETA_PHI),
		"Value of beta_phi was not within tolerance. Got: " << alb_test.beta_phi << " Exp: " << BETA_PHI);
	BOOST_CHECK_MESSAGE(test_val(alb_test.alpha_phi, ALPHA_PHI) && test_val(alb_test.beta_phi, BETA_PHI),
		"Value of theta objective func at estimated (theta) hparams: " << bobquant.eval_theta_objective_func(alb_test) <<
		"\nAt expected hparams: " << bobquant.eval_theta_objective_func(expected) <<
		"\nValue of phi objective func at estimated hparams: " << bobquant.eval_phi_objective_func(alb_test) <<
		"\nAt expected hparams: " << bobquant.eval_phi_objective_func(expected));
	/*BOOST_CHECK_MESSAGE(test_val(alb_test.alpha/alb_test.beta, ALPHA_PHI/BETA_PHI),
		"Ratio of alpha and beta was not within tolerance: Got: " << alb_test.alpha/alb_test.beta <<
		"(a=" << alb_test.alpha << ", b=" << alb_test.beta << ") Exp: " << ALPHA_PHI/BETA_PHI <<
		"(a=" << ALPHA_PHI << ", b=" << BETA_PHI << ")");*/
}

BOOST_AUTO_TEST_SUITE(test)
//
// component function tests
//

BOOST_AUTO_TEST_CASE(static_functions) {
	const long int XIJ = 1;
	const long int NIJ = 6;
	const long int XCJ = 3;
	const long int NCJ = 7;
	const double ALPHA_THETA = 1;
	const double BETA_THETA = 10;
	const double ALPHA_PHI = 1;
	const double BETA_PHI = 100;
	const double LP_XIJ_THETA_STD = -3.688879;
	const double LP_XCJ_PHI_STD = -8.026631;
	const double DLP_XIJ_DALPHA_STD = -2.380729;
	const double DLP_XIJ_DBETA_STD = -0.1291667;
	double lp_xij_theta = hts::bayes_orient_bias_quant_f::theta_hparams_optimizable::lp_xij_given_hparams(
		XIJ, NIJ, ALPHA_THETA, BETA_THETA);
	double lp_xcj_phi = hts::bayes_orient_bias_quant_f::phi_hparams_optimizable::lp_xcj_given_hparams(
		XCJ, NCJ, ALPHA_THETA, BETA_THETA, ALPHA_PHI, BETA_PHI);
	double dlp_xij_dalpha = hts::bayes_orient_bias_quant_f::theta_hparams_optimizable::dlp_xij_given_hparams_dalpha(
		XIJ, NIJ, ALPHA_THETA, BETA_THETA);
	double dlp_xij_dbeta = hts::bayes_orient_bias_quant_f::theta_hparams_optimizable::dlp_xij_given_hparams_dbeta(
		XIJ, NIJ, ALPHA_THETA, BETA_THETA);
	BOOST_CHECK_MESSAGE(test_val(lp_xij_theta, LP_XIJ_THETA_STD),
		"lp_xij_given_hparams (theta) not within tolerance: Got: " << lp_xij_theta << " Exp: " << LP_XIJ_THETA_STD);
	BOOST_CHECK_MESSAGE(test_val(lp_xcj_phi, LP_XCJ_PHI_STD),
		"lp_xcj_given_hparams (phi) not within tolerance: Got: " << lp_xcj_phi << " Exp: " << LP_XCJ_PHI_STD);
	BOOST_CHECK_MESSAGE(test_val(dlp_xij_dalpha, DLP_XIJ_DALPHA_STD),
		"dlp_xij_given_hparams_dalpha (theta) not within tolerance: Got: " << dlp_xij_dalpha << " Exp: " << DLP_XIJ_DALPHA_STD);
	BOOST_CHECK_MESSAGE(test_val(dlp_xij_dbeta, DLP_XIJ_DBETA_STD),
		"dlp_xij_given_hparams_dbeta (theta) not within tolerance: Got: " << dlp_xij_dbeta << " Exp: " << DLP_XIJ_DBETA_STD);
}

BOOST_AUTO_TEST_CASE(member_functions) {
	const double THETA_OBJFUNC_STD = -21.26166;
	const double PHI_OBJFUNC_STD = -31.12604;
	const double THETA_GRAD_DALPHA_STD = -0.2056439;
	const double THETA_GRAD_DBETA_STD = 0.03687712;
	hts::hparams TEST;
	TEST.alpha_theta = 1;
	TEST.beta_theta = 10;
	TEST.alpha_phi = 1;
	TEST.beta_phi = 100;
	const long int XC_ARR[] = {2, 0, 0, 2, 1, 2, 3, 0, 1, 0, 0, 0, 1, 2, 1, 0, 1, 1, 4, 2};
	const long int XI_ARR[] = {0, 0, 0, 1, 2, 0, 1, 2, 1, 0, 0, 2, 0, 0, 0, 2, 0, 0, 1, 0};
	const long int NC_ARR[] = {9, 9, 10, 6, 10, 5, 7, 7, 7, 10, 9, 7, 6, 6, 7, 7, 7, 7, 7, 10};
	const long int NI_ARR[] = {8, 9, 9, 6, 9, 5, 6, 7, 6, 9, 8, 6, 5, 5, 7, 6, 7, 7, 7, 9};
	vector<long int> xc_vec (XC_ARR, XC_ARR + sizeof(XC_ARR)/sizeof(long int));
	vector<long int> xi_vec (XI_ARR, XI_ARR + sizeof(XI_ARR)/sizeof(long int));
	vector<long int> nc_vec (NC_ARR, NC_ARR + sizeof(NC_ARR)/sizeof(long int));
	vector<long int> ni_vec (NI_ARR, NI_ARR + sizeof(NI_ARR)/sizeof(long int));
	hts::bayes_orient_bias_quant_f bobquant (nuc_T, nuc_G);
	bobquant.set_data(xc_vec, xi_vec, nc_vec, ni_vec);
	// these are negated because the objective functions are negated internally
	double theta_objfunc = -bobquant.eval_theta_objective_func(TEST);
	double phi_objfunc = -bobquant.eval_phi_objective_func(TEST);
	double theta_grad_dalpha = bobquant.eval_theta_grad_dalpha(TEST);
	double theta_grad_dbeta = bobquant.eval_theta_grad_dbeta(TEST);
	BOOST_CHECK_MESSAGE(test_val(theta_objfunc, THETA_OBJFUNC_STD),
		"lp_xi_given_hparams (theta objective function) not within tolerance. Got: " << theta_objfunc << " Exp: " << THETA_OBJFUNC_STD);
	BOOST_CHECK_MESSAGE(test_val(theta_grad_dalpha, THETA_GRAD_DALPHA_STD),
		"dlp_xi_given_hparams_dalpha (theta gradient wrt alpha) not within tolerance. Got: " << theta_grad_dalpha << " Exp: " << THETA_GRAD_DALPHA_STD);
	BOOST_CHECK_MESSAGE(test_val(theta_grad_dbeta, THETA_GRAD_DBETA_STD),
		"dlp_xi_given_hparams_dbeta (theta gradient wrt beta) not within tolerance. Got: " << theta_grad_dbeta << " Exp: " << THETA_GRAD_DBETA_STD);
	BOOST_CHECK_MESSAGE(test_val(phi_objfunc, PHI_OBJFUNC_STD),
		"lp_xc_given_hparams (phi objective function) not within tolerance. Got: " << phi_objfunc << " Exp: " << PHI_OBJFUNC_STD);
}

//
// freq tests
//
BOOST_AUTO_TEST_CASE(freq_hs_hd) {
	const size_t N = FREQ_N_READS;
	const double THETA = .1;
	const double PHI = .1;
	freq_obquant_sim_test(N, THETA, PHI);
}

BOOST_AUTO_TEST_CASE(freq_hs_ld) {
	const size_t N = FREQ_N_READS;
	const double THETA = .1;
	const double PHI = .01;
	freq_obquant_sim_test(N, THETA, PHI);
}

BOOST_AUTO_TEST_CASE(freq_ls_hd) {
	const size_t N = FREQ_N_READS;
	const double THETA = .01;
	const double PHI = .1;
	freq_obquant_sim_test(N, THETA, PHI);
}

BOOST_AUTO_TEST_CASE(freq_ms_hd) {
	const size_t N = FREQ_N_READS;
	const double THETA = .05;
	const double PHI = .1;
	freq_obquant_sim_test(N, THETA, PHI);
}

BOOST_AUTO_TEST_CASE(freq_ms_ld) {
	const size_t N = FREQ_N_READS;
	const double THETA = .05;
	const double PHI = .01;
	freq_obquant_sim_test(N, THETA, PHI);
}

BOOST_AUTO_TEST_CASE(freq_ls_ld) {
	const size_t N = FREQ_N_READS;
	const double THETA = .01;
	const double PHI = .01;
	freq_obquant_sim_test(N, THETA, PHI);
}

//
// bayes tests
//
BOOST_AUTO_TEST_CASE(bayes_hs_hd) {
	const double ALPHA_THETA = 1;
	const double BETA_THETA = 10;
	const double ALPHA_PHI = 1;
	const double BETA_PHI = 10;
	bayes_obquant_sim_test(BAYES_N_LOCS, BAYES_MIN_RCOUNT, BAYES_MAX_RCOUNT, ALPHA_THETA, BETA_THETA, ALPHA_PHI, BETA_PHI);
}

BOOST_AUTO_TEST_CASE(bayes_hs_ld) {
	const double ALPHA_THETA = 1;
	const double BETA_THETA = 10;
	const double ALPHA_PHI = 1;
	const double BETA_PHI = 100;
	bayes_obquant_sim_test(BAYES_N_LOCS, BAYES_MIN_RCOUNT, BAYES_MAX_RCOUNT, ALPHA_THETA, BETA_THETA, ALPHA_PHI, BETA_PHI);
}

BOOST_AUTO_TEST_CASE(bayes_ls_hd) {
	const double ALPHA_THETA = 1;
	const double BETA_THETA = 100;
	const double ALPHA_PHI = 1;
	const double BETA_PHI = 10;
	bayes_obquant_sim_test(BAYES_N_LOCS, BAYES_MIN_RCOUNT, BAYES_MAX_RCOUNT, ALPHA_THETA, BETA_THETA, ALPHA_PHI, BETA_PHI);
}

BOOST_AUTO_TEST_CASE(bayes_ms_hd) {
	const double ALPHA_THETA = 1;
	const double BETA_THETA = 50;
	const double ALPHA_PHI = 1;
	const double BETA_PHI = 10;
	bayes_obquant_sim_test(BAYES_N_LOCS, BAYES_MIN_RCOUNT, BAYES_MAX_RCOUNT, ALPHA_THETA, BETA_THETA, ALPHA_PHI, BETA_PHI);
}

BOOST_AUTO_TEST_CASE(bayes_ms_ld) {
	const double ALPHA_THETA = 1;
	const double BETA_THETA = 50;
	const double ALPHA_PHI = 1;
	const double BETA_PHI = 100;
	bayes_obquant_sim_test(BAYES_N_LOCS, BAYES_MIN_RCOUNT, BAYES_MAX_RCOUNT, ALPHA_THETA, BETA_THETA, ALPHA_PHI, BETA_PHI);
}

BOOST_AUTO_TEST_CASE(bayes_1s_1d) {
	const double ALPHA_THETA = 1;
	const double BETA_THETA = 100;
	const double ALPHA_PHI = 1;
	const double BETA_PHI = 100;
	bayes_obquant_sim_test(BAYES_N_LOCS, BAYES_MIN_RCOUNT, BAYES_MAX_RCOUNT, ALPHA_THETA, BETA_THETA, ALPHA_PHI, BETA_PHI);
}

BOOST_AUTO_TEST_SUITE_END()