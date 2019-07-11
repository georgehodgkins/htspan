#define BOOST_TEST_MODULE orient_bias_quant_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstddef>

#include <alglib/alglibmisc.h>

#include "../htspan/nucleotide.hpp"
#include "../htspan/orient_bias_quant.hpp"

#define TEST_EPS .01
#define TEST_RAT .2

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
	hts::bayes_orient_bias_quant_f bobquant (nuc_T, nuc_C);
	bobquant.simulate(ns_vec, ALPHA_THETA, BETA_THETA, ALPHA_PHI, BETA_PHI);
	hts::alpha_beta alb_test = bobquant(BAYES_BSIZE, BAYES_NEPOCHS, BAYES_LEARN_RATE, BAYES_EPS,
		BAYES_INIT_ALPHA, BAYES_INIT_BETA, BAYES_INIT_ALPHA, BAYES_INIT_BETA);
	BOOST_CHECK_MESSAGE(test_val(alb_test.alpha/alb_test.beta, ALPHA_PHI/BETA_PHI),
		"Ratio of alpha and beta was not within tolerance: Got: " << alb_test.alpha/alb_test.beta <<
		"(a=" << alb_test.alpha << ", b=" << alb_test.beta << ") Exp: " << ALPHA_PHI/BETA_PHI <<
		"(a=" << ALPHA_PHI << ", b=" << BETA_PHI << ")");
}

BOOST_AUTO_TEST_SUITE(test)

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