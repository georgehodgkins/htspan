#define BOOST_TEST_MODULE orient_bias_quant_module
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#define TEST_EPS 1e-3

#include <cmath>

#include "../htspan/nucleotide.hpp"
#include "../htspan/orient_bias_quant.hpp"

void common_obquant_operator_test(const int SET_XC, const int SET_NC, const int SET_XI, const int SET_NI,
		const double PHI_STD) {
	hts::orient_bias_quant_f obquant (nuc_G, nuc_T);
	obquant.xc = SET_XC;
	obquant.nc = SET_NC;
	obquant.xi = SET_XI;
	obquant.ni = SET_NI;
	double result = obquant();
	BOOST_CHECK_MESSAGE(abs(result - PHI_STD) < TEST_EPS,
		"Phi estimate for quantification: got: " << result << ", expected: " << PHI_STD);
}

BOOST_AUTO_TEST_SUITE (orient_bias_quant)

BOOST_AUTO_TEST_CASE (operator_only_1) {
	const int SET_XC = 5;
	const int SET_NC = 100;
	const int SET_XI = 2;
	const int SET_NI = 90;
	const double PHI_STD = 0.02840909;
	common_obquant_operator_test(SET_XC, SET_NC, SET_XI, SET_NI, PHI_STD);
}

BOOST_AUTO_TEST_CASE (operator_only_2) {
	const int SET_XC = 1;
	const int SET_NC = 100;
	const int SET_XI = 1;
	const int SET_NI = 90;
	const double PHI_STD = 0.0;
	common_obquant_operator_test(SET_XC, SET_NC, SET_XI, SET_NI, PHI_STD);
}

BOOST_AUTO_TEST_SUITE_END();
