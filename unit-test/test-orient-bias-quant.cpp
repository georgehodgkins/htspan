#define BOOST_TEST_MODULE orient_bias_quant_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdlib>
#include <fstream>

#include <alglib/alglibmisc.h>

#include "htspan/nucleotide.hpp"
#include "htspan/freq_orient_bias_quant.hpp"
#include "htspan/bayes_orient_bias_quant.hpp"

// set slightly looser standards than the default
#define TEST_RAT .1
#define TEST_EPS 1e-6
#include "test.hpp"

#define FREQ_N_READS 75000

#define BAYES_STEPPER stograd::stepper::adam<double>
#define BAYES_LEARN_RATE 1e-4
#define BAYES_EPS 1e-6  // avoid early convergence to local minimum
#define BAYES_BSIZE 100
#define BAYES_NEPOCHS 3000
#define BAYES_INIT_ALPHA 1
#define BAYES_INIT_BETA 1
#define BAYES_N_LOCS 75000
#define BAYES_MIN_RCOUNT 50
#define BAYES_MAX_RCOUNT 150

using namespace std;

// Boost fixture to allow access to command line arguments
struct seed_fixture {
	int SIM_SEED;

	seed_fixture()
		{
			if (boost::unit_test::framework::master_test_suite().argc > 1) {
				SIM_SEED = atoi(boost::unit_test::framework::master_test_suite().argv[1]);
			} else {
				SIM_SEED = 0;
			}
		}

	~seed_fixture() {}
};

void freq_obquant_sim_test (const size_t N, const double THETA, const double PHI) {
	hts::freq_orient_bias_quant_f fobquant (nuc_T, nuc_C);
	fobquant.simulate(N, THETA, PHI);
	double phi_test = fobquant();
	BOOST_CHECK_MESSAGE(test_val(phi_test, PHI),
		"Phi value returned by estimator was not within tolerance. Got: " << phi_test << " Exp: " << PHI);
}

void bayes_obquant_sim_test (const size_t N, const size_t LOC_RD_MIN, const size_t LOC_RD_MAX,
		const double ALPHA_THETA_STD, const double BETA_THETA_STD, const double ALPHA_PHI_STD, const double BETA_PHI_STD, int sim_seed) { 
	// generate random read counts within the requested range
	if (sim_seed == 0) {	
		alglib::hqrndstate rng;
		alglib::hqrndrandomize(rng);
		sim_seed = alglib::hqrnduniformi(rng, 1 << 16); // uniform int between 0 and 2^15
	}
	BOOST_TEST_MESSAGE("Bayesian quant model test with true hparams alpha_theta = " << ALPHA_THETA_STD <<
		", beta_theta = " << BETA_THETA_STD << ", alpha_phi = " << ALPHA_PHI_STD << ", beta_phi = " << BETA_PHI_STD <<
		"\nSimulating " << N << " loci with read counts between " << LOC_RD_MIN << " and " << LOC_RD_MAX <<  " (seed = " << sim_seed << ")");
	// simulate oberved variables
	double THETA_MEAN = ALPHA_THETA_STD/(ALPHA_THETA_STD + BETA_THETA_STD);
	double PHI_MEAN = ALPHA_PHI_STD/(ALPHA_PHI_STD + BETA_PHI_STD);
	hts::bayes_orient_bias_quant_f bobquant (nuc_T, nuc_C);
	bobquant.simulate(N, ALPHA_THETA_STD, BETA_THETA_STD, ALPHA_PHI_STD, BETA_PHI_STD, LOC_RD_MIN, LOC_RD_MAX, sim_seed);
	hts::freq_orient_bias_quant_f fobquant (nuc_T, nuc_C);
	fobquant.copy_data(bobquant.m.xc_vec, bobquant.m.xi_vec, bobquant.m.nc_vec, bobquant.m.ni_vec);
	BOOST_CHECK_MESSAGE(test_val(fobquant.theta_hat(), THETA_MEAN),
		"Frequentist estimate of generated theta: Got: " << fobquant.theta_hat() << " Exp: " << THETA_MEAN);
	BOOST_CHECK_MESSAGE(test_val(fobquant(), PHI_MEAN),
		"Frequentist estimate of generated phi: Got: " << fobquant() << " Exp: " << PHI_MEAN);
	// run the test
	hts::hparams got = bobquant.operator()<BAYES_STEPPER>(BAYES_BSIZE, BAYES_NEPOCHS, BAYES_LEARN_RATE, BAYES_EPS,
		BAYES_INIT_ALPHA, BAYES_INIT_BETA, BAYES_INIT_ALPHA, BAYES_INIT_BETA);
	double E_theta = got.alpha_theta/(got.alpha_theta + got.beta_theta);
	double E_phi = got.alpha_phi/(got.alpha_phi + got.beta_phi);
	BOOST_TEST_MESSAGE("Theta estimate epoch count: " << bobquant.n_theta_epochs);
	BOOST_CHECK_MESSAGE(test_val(got.alpha_theta, ALPHA_THETA_STD), 
		"alpha_theta: Got: " << got.alpha_theta << " (log = " << log(got.alpha_theta) <<
		") Exp: " << ALPHA_THETA_STD << " (log = " << log(ALPHA_THETA_STD) << ")");
	BOOST_CHECK_MESSAGE(test_val(got.beta_theta, BETA_THETA_STD),
		"beta_theta: Got: " << got.beta_theta << " (log = " << log(got.beta_theta) <<
		") Exp: " << BETA_THETA_STD << " (log = " << log(BETA_THETA_STD) << ")");
	BOOST_CHECK_MESSAGE(test_val(E_theta, THETA_MEAN),
		"E[theta]: Got: " << E_theta << " Exp: " << THETA_MEAN);
	BOOST_TEST_MESSAGE("Phi estimate epoch count: " << bobquant.n_phi_epochs);
	BOOST_CHECK_MESSAGE(test_val(got.alpha_phi, ALPHA_PHI_STD), 
		"alpha_phi: Got: " << got.alpha_phi << " (log = " << log(got.alpha_phi) <<
		") Exp: " << ALPHA_PHI_STD << " (log = " << log(ALPHA_PHI_STD) << ")");
	BOOST_CHECK_MESSAGE(test_val(got.beta_phi, BETA_PHI_STD),
		"beta_phi: Got: " << got.beta_phi << " (log = " << log(got.beta_phi) <<
		") Exp: " << BETA_PHI_STD << " (log = " << log(BETA_PHI_STD) << ")");
	BOOST_CHECK_MESSAGE(test_val(E_phi, PHI_MEAN),
		"E[phi]: Got: " << E_phi << " Exp: " << PHI_MEAN);
}

BOOST_FIXTURE_TEST_SUITE(test, seed_fixture)
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
	const double LP_XCJ_PHI_STD = -3.198318;
	const double DLP_XIJ_DALPHA_STD = 0.5482393;
	const double DLP_XIJ_DBETA_STD = -0.02916667;
	double lp_xcj_phi = hts::bayes_quant_model::lp_xcj_given_hparams(
		XCJ, NCJ, ALPHA_THETA, BETA_THETA, ALPHA_PHI, BETA_PHI);
	double dlp_xij_dalpha = hts::bayes_quant_model::dlp_xij_given_hparams_dalpha(
		XIJ, NIJ, ALPHA_THETA, BETA_THETA);
	double dlp_xij_dbeta = hts::bayes_quant_model::dlp_xij_given_hparams_dbeta(
		XIJ, NIJ, ALPHA_THETA, BETA_THETA);
	BOOST_CHECK_MESSAGE(test_val(lp_xcj_phi, LP_XCJ_PHI_STD),
		"lp_xcj_given_hparams (phi) not within tolerance: Got: " << lp_xcj_phi << " Exp: " << LP_XCJ_PHI_STD);
	BOOST_CHECK_MESSAGE(test_val(dlp_xij_dalpha, DLP_XIJ_DALPHA_STD),
		"dlp_xij_given_hparams_dalpha (theta) not within tolerance: Got: " << dlp_xij_dalpha << " Exp: " << DLP_XIJ_DALPHA_STD);
	BOOST_CHECK_MESSAGE(test_val(dlp_xij_dbeta, DLP_XIJ_DBETA_STD),
		"dlp_xij_given_hparams_dbeta (theta) not within tolerance: Got: " << dlp_xij_dbeta << " Exp: " << DLP_XIJ_DBETA_STD);
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
	const double ALPHA_THETA_STD = 1;
	const double BETA_THETA_STD = 10;
	const double ALPHA_PHI_STD = 1;
	const double BETA_PHI_STD = 10;
	bayes_obquant_sim_test(BAYES_N_LOCS, BAYES_MIN_RCOUNT, BAYES_MAX_RCOUNT, ALPHA_THETA_STD, BETA_THETA_STD, ALPHA_PHI_STD, BETA_PHI_STD, SIM_SEED);
}

BOOST_AUTO_TEST_CASE(bayes_hs_ld) {
	const double ALPHA_THETA_STD = 1;
	const double BETA_THETA_STD = 10;
	const double ALPHA_PHI_STD = 1;
	const double BETA_PHI_STD = 100;
	bayes_obquant_sim_test(BAYES_N_LOCS, BAYES_MIN_RCOUNT, BAYES_MAX_RCOUNT, ALPHA_THETA_STD, BETA_THETA_STD, ALPHA_PHI_STD, BETA_PHI_STD, SIM_SEED);
}

BOOST_AUTO_TEST_CASE(bayes_ls_hd) {
	const double ALPHA_THETA_STD = 1;
	const double BETA_THETA_STD = 100;
	const double ALPHA_PHI_STD = 1;
	const double BETA_PHI_STD = 10;
	bayes_obquant_sim_test(BAYES_N_LOCS, BAYES_MIN_RCOUNT, BAYES_MAX_RCOUNT, ALPHA_THETA_STD, BETA_THETA_STD, ALPHA_PHI_STD, BETA_PHI_STD, SIM_SEED);
}

BOOST_AUTO_TEST_CASE(bayes_ms_hd) {
	const double ALPHA_THETA_STD = 1;
	const double BETA_THETA_STD = 50;
	const double ALPHA_PHI_STD = 1;
	const double BETA_PHI_STD = 10;
	bayes_obquant_sim_test(BAYES_N_LOCS, BAYES_MIN_RCOUNT, BAYES_MAX_RCOUNT, ALPHA_THETA_STD, BETA_THETA_STD, ALPHA_PHI_STD, BETA_PHI_STD, SIM_SEED);
}

BOOST_AUTO_TEST_CASE(bayes_ms_ld) {
	const double ALPHA_THETA_STD = 1;
	const double BETA_THETA_STD = 50;
	const double ALPHA_PHI_STD = 1;
	const double BETA_PHI_STD = 100;
	bayes_obquant_sim_test(BAYES_N_LOCS, BAYES_MIN_RCOUNT, BAYES_MAX_RCOUNT, ALPHA_THETA_STD, BETA_THETA_STD, ALPHA_PHI_STD, BETA_PHI_STD, SIM_SEED);
}

BOOST_AUTO_TEST_CASE(bayes_ls_ld) {
	const double ALPHA_THETA_STD = 1;
	const double BETA_THETA_STD = 100;
	const double ALPHA_PHI_STD = 1;
	const double BETA_PHI_STD = 100;
	bayes_obquant_sim_test(BAYES_N_LOCS, BAYES_MIN_RCOUNT, BAYES_MAX_RCOUNT, ALPHA_THETA_STD, BETA_THETA_STD, ALPHA_PHI_STD, BETA_PHI_STD, SIM_SEED);
}

BOOST_AUTO_TEST_SUITE_END()
