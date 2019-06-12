#define BOOST_TEST_MODULE bayes_orient_bias_filter_module
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#define TEST_EPS 1e-3

#include <cmath>

#include <vector>

#include "htspan/bayes_orient_bias_filter.hpp"

using namespace std;

// UNDER CONSTRUCTION
// TODO: add specific test cases

void common_cgrid_test(hts::bayes_orient_bias_filter_f& bobfilter, const double GRID_PARAM, const double GRID_EPS,
		const double GRID_4F_STD, const double GRID_4L_STD, const size_t GRID_SIZE_STD) {
	BOOST_TEST_MESSAGE("Running centered grid test with parameter " << GRID_PARAM << ":");
	vector<double> test_grid = bobfilter.generate_cgrid(GRID_PARAM, GRID_EPS);
	BOOST_CHECK_MESSAGE(test_grid.front() == GRID_EPS,
		"First member is not equal to set eps (got: " << test_grid.front() << ")");
	BOOST_CHECK_MESSAGE(test_grid.back() == (1 - GRID_EPS),
		"Last member is not equal to 1 - set eps (got: " << test_grid.back() << ")");
	BOOST_CHECK_MESSAGE(test_grid.size() == GRID_SIZE_STD,
		"Grid size: got: " << test_grid.size() << ", expected: " << GRID_SIZE_STD);
	BOOST_CHECK_MESSAGE(abs(test_grid[3] - GRID_4F_STD) < TEST_EPS, 
		"Grid[3] (4th member): got: " << test_grid[3] << ", expected: " << GRID_4F_STD);
	BOOST_CHECK_MESSAGE(abs(test_grid[grid.size() - 4] - GRID_4L_STD) < TEST_EPS,
		"Grid[-4] (4th from last member): got: " << test_grid[grid.size() - 4] << ", expected:" << GRID_4L_STD);
}

void common_integration_test(const char TSVNAME[], const double PHI_0, const double THETA_0,
		const double ALPHA_PHI, const double BETA_PHI, const double PRIOR_ALT,
		const double PHI_INT_STD, const double THETA_INT_STD, const double EV_NULL_STD, 
		const double EV_ALT_STD, const double LPOSTERIOR_STD) {
	BOOST_TEST_CHECKPOINT("Initializing filter and reading in data");
	hts::bayes_orient_bias_filter bobfilter(nuc_G, nuc_T, 200);
	bobfilter.read(TSVNAME);
	BOOST_TEST_CHECKPOINT("Evaluating phi integrand");
	bobfilter.theta_t = THETA_0;
	bobfilter.alpha_phi = ALPHA_PHI;
	bobfilter.beta_phi = BETA_PHI;
	double phi_int = bobfilter.phi_integrand(PHI_0);
	BOOST_CHECK_MESSAGE(abs(phi_int - PHI_INT_STD) < TEST_EPS,
		"Phi integrand: got: " << phi_int << ", expected: " << PHI_INT_STD);
	BOOST_TEST_CHECKPOINT("Evaluating theta integrand");
	double theta_int = bobfilter.theta_integrand(THETA_0);
	BOOST_CHECK_MESSAGE(abs(theta_int - THETA_INT_STD) < TEST_EPS,
		"Theta integrand: got: " << theta_int << ", expected: " << THETA_INT_STD);
	BOOST_TEST_CHECKPOINT("Evaluating evidence");
	hts::evidence_rtn ev = bobfilter.model_evidence(ALPHA_PHI, BETA_PHI);
	BOOST_CHECK_MESSAGE(abs(ev.null - EV_NULL_STD) < TEST_EPS,
		"Evidence for null model: got: " << ev.null << ", expected: " << EV_NULL_STD);
	BOOST_CHECK_MESSAGE(abs(ev.alt - EV_ALT_STD) < TEST_EPS,
		"Evidence for alternate model: got: " << ev.alt << ", expected: " << EV_ALT_STD);
	BOOST_TEST_CHECKPOINT("Evaluating posterior probability");
	double lposterior = bobfilter(PRIOR_ALT, ALPHA_PHI, BETA_PHI);
	BOOST_CHECK_MESSAGE(abs(lposterior - LPOSTERIOR_STD) < TEST_EPS,
		"Log posterior probability: got: " << lposterior << ", expected: " << LPOSTERIOR_STD);
} 
