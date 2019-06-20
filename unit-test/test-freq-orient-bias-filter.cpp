#define BOOST_TEST_MODULE freq_orient_bias_filter_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cmath>

#include "htspan/freq_orient_bias_filter.hpp"
#include "htspan/io/snv.hpp"
#include "htspan/nucleotide.hpp"

#include "test.hpp"

void common_filter_math_test (const char TSVNAME[], const double THETA_0, const double PHI_0,
		const double THETA_F_STD, const double PVAL_F_STD,
		const double THETA_CA_STD, const double PHI_CA_STD, const double PVAL_CA_STD) {
	// basic initialization
	BOOST_TEST_CHECKPOINT("Initializing data struct");
	hts::orient_bias_data data(nuc_G, nuc_T, 200);
	BOOST_TEST_CHECKPOINT("Reading data into data struct");
	data.read(TSVNAME);
	BOOST_TEST_CHECKPOINT("Estimating theta for fixed phi case");
	hts::freq_orient_bias_filter_f fobfilter (data);
	// fixed phi cases
	double theta_fixed_phi = fobfilter.estimate_theta_given(PHI_0, THETA_0);
	// pass this test if minimizer outperforms the R minimizer
	BOOST_CHECK_MESSAGE(test_val(theta_fixed_phi, THETA_F_STD) || 
		(-fobfilter.lp_bases_given(theta_fixed_phi, PHI_0) < -fobfilter.lp_bases_given(THETA_F_STD, PHI_0)), 
		"Theta estimate for fixed phi: got:" << theta_fixed_phi << ", expected: " << THETA_F_STD << "\n"
		<< "Objective f(got_theta, phi_0): " << -fobfilter.lp_bases_given(theta_fixed_phi, PHI_0) <<
		", f(exp_theta, phi_0): " << -fobfilter.lp_bases_given(THETA_F_STD, PHI_0));
	BOOST_TEST_CHECKPOINT("Calculating p-val for fixed phi case");
	double pval_fixed_phi = fobfilter(PHI_0);
	BOOST_CHECK_MESSAGE(test_val(pval_fixed_phi, PVAL_F_STD), 
		"P-val for fixed phi: got:" << pval_fixed_phi << ", expected: " << PVAL_F_STD);
	// unknown phi cases
	BOOST_TEST_CHECKPOINT("Estimating theta and phi with unknown phi");
	hts::theta_and_phi unk_phi = fobfilter.estimate_theta_phi(THETA_0, PHI_0);
	BOOST_CHECK_MESSAGE(test_val(unk_phi.theta, THETA_CA_STD),
		"Theta estimate with unknown phi: got:" << unk_phi.theta << ", expected: " << THETA_CA_STD);
	BOOST_CHECK_MESSAGE(test_val(unk_phi.phi, PHI_CA_STD),
		"Phi estimate with unknown phi: got:" << unk_phi.phi << ", expected: " << PHI_CA_STD);
	// pass the test if the minimizer outperforms the R minimizer
	BOOST_CHECK_MESSAGE((test_val(unk_phi.theta, THETA_CA_STD) && test_val(unk_phi.phi, PHI_CA_STD)) || 
		-fobfilter.lp_bases_given(unk_phi.theta, unk_phi.phi) < -fobfilter.lp_bases_given(THETA_CA_STD, PHI_CA_STD),
		"For unknown phi: Objective f(got_theta, got_phi): " << -fobfilter.lp_bases_given(unk_phi.theta, unk_phi.phi) <<
		", f(exp_theta, exp_phi): " << -fobfilter.lp_bases_given(THETA_CA_STD, PHI_CA_STD));
	BOOST_TEST_CHECKPOINT("Calculating p-val with unknown phi estimation");
	double pval_coord_ascent = fobfilter(PHI_0, false);
	BOOST_CHECK_MESSAGE(test_val(pval_coord_ascent, PVAL_CA_STD), 
		"P-val for unknown phi: got:" << pval_coord_ascent << ", expected: " << PVAL_CA_STD);
}


BOOST_AUTO_TEST_SUITE(test)

BOOST_AUTO_TEST_CASE (stat_reader) {
	const char TSVNAME[] = "../sim-data/obrs_low-signal_high-damage_data.tsv";
	const size_t DLEN = 200;
	const double ERR_0 = 0.000781619442268712;
	const int8_t BASE_0 = 0;
	const bool ORIENT_0 = true;
	const double ERR_99 = 0.000691675429225815;
	const int8_t BASE_99 = 1;
	const bool ORIENT_99 = true;
	const double ERR_199 = 0.00132148527516072;
	const int8_t BASE_199 = 1;
	const bool ORIENT_199 = true;
	BOOST_TEST_MESSAGE("Running stat reader test:");
	hts::orient_bias_data data(nuc_T, nuc_G, 0);
	data.read(TSVNAME);
	BOOST_CHECK_MESSAGE((data.bases.size() == DLEN && data.orients.size() == DLEN && data.errors.size() == DLEN),
		"Stat reader did not read the correct number of reads. Read: " << data.bases.size());
	BOOST_CHECK_MESSAGE(abs(data.errors[0] - ERR_0) < TEST_EPS,
		"0th member of errors vector does not match.");
	BOOST_CHECK_MESSAGE(data.bases[99] == BASE_99,
		"99th member of bases vector does not match.");
	BOOST_CHECK_MESSAGE(abs(data.errors[99] - ERR_99) < TEST_EPS,
		"99th member of errors vector does not match.");
	BOOST_CHECK_MESSAGE(data.orients[199] == ORIENT_199,
		"199th member of orients vector does not match.");
	BOOST_CHECK_MESSAGE(abs(data.errors[199] - ERR_199) < TEST_EPS,
		"199th member of errors vector does not match.");
}

BOOST_AUTO_TEST_CASE (tsv_reader) {
	const char SNVNAME[] = "../sim-data/../../data/snv.tsv";
	const size_t SLEN = 3;
	const nuc_t REF_F = nuc_C;
	const long int POS_F = 7674420;
	const nuc_t REF_L = nuc_C;
	const long int POS_L = 7674361;
	BOOST_TEST_MESSAGE("Running TSV reader test: ");
	hts::snv::tsv_reader snvr (SNVNAME, NULL);
	hts::snv::record recF, recL;
	snvr.next(recF);
	while (snvr.next(recL)) {}
	BOOST_CHECK_MESSAGE(recF.nt_ref == REF_F,
		"First record ref nucleotide does not match.");
	BOOST_CHECK_MESSAGE(recF.pos == POS_F-1,
		"First record position does not match. Got: " << recF.pos << ", expected: " << POS_F-1);
	BOOST_CHECK_MESSAGE(recL.nt_ref == REF_L,
		"Last record ref nucleotide does not match.");
	BOOST_CHECK_MESSAGE(recL.pos == POS_L-1,
		"Last record position does not match. Got: " << recL.pos << ", expected: " << POS_L-1);
}

BOOST_AUTO_TEST_CASE (lp_functions) {
	const int8_t BASES[] = {0, 1, 2, 0, 0, 1, 2, 0, 1, 1};
	const double ERRS[] = {1e-2, 1e-3, 1e-4, 1e-3, 1e-3, 1e-2, 1e-5, 1e-2, 1e-3, 1e-4};
	const bool ORIENTS[] = {true, true, false, true, true, false, true, false, true, false};
	const double PHI = .04;
	const double THETA = .01;
	const double LP_REF_ERR = 1e-5;
	const double LP_REF_STD = -0.05046151;
	const double LP_ALT_ERR = 1e-5;
	const double LP_ALT_STD = -3.011807;
	const double LP_OTH_ERR = 1e-4;
	const double LP_OTH_STD = -3.218176;
	const double LP_FULL_STD = -28.65358;
	const double LP_BASE_STD = -36.10886;
	const double DEV_STD = 14.91057;
	BOOST_TEST_MESSAGE("Running log-prob function tests:");
	std::vector<int8_t> base_t (BASES, BASES + sizeof(BASES)/sizeof(int8_t));
	std::vector<double> err_t (ERRS, ERRS + sizeof(ERRS)/sizeof(double));
	std::vector<bool> orient_t (ORIENTS, ORIENTS + sizeof(ORIENTS)/sizeof(bool));
	hts::orient_bias_data data (nuc_G, nuc_T, 0);
	data.bases.swap(base_t);
	data.errors.swap(err_t);
	data.orients.swap(orient_t);
	hts::freq_orient_bias_filter_f fobfilter(data);
	double lp_ref = fobfilter.lp_ref_base_given(THETA, PHI, LP_REF_ERR);
	BOOST_CHECK_MESSAGE((abs(lp_ref - LP_REF_STD) < TEST_EPS),
		"LP with ref base: got: " << lp_ref << ", expected: " << LP_REF_STD);
	double lp_alt = fobfilter.lp_alt_base_given(THETA, PHI, LP_ALT_ERR);
	BOOST_CHECK_MESSAGE((abs(lp_alt - LP_ALT_STD) < TEST_EPS),
		"LP with alt base: got: " << lp_alt << ", expected: " << LP_ALT_STD);
	double lp_oth = fobfilter.lp_other_base_given(THETA, PHI, LP_OTH_ERR);
	BOOST_CHECK_MESSAGE((abs(lp_oth - LP_OTH_STD) < TEST_EPS),
		"LP with oth base: got: " << lp_alt << ", expected: " << LP_OTH_STD);
	double lp_full = fobfilter.lp_bases_given(THETA, PHI);
	BOOST_CHECK_MESSAGE((abs(lp_full - LP_FULL_STD) < TEST_EPS),
		"LP with test data: got: " << lp_full << ", expected: " << LP_FULL_STD);
	double lp_base = fobfilter.lp_bases_given(0, PHI);
	BOOST_CHECK_MESSAGE((abs(lp_base - LP_BASE_STD) < TEST_EPS),
		"LP for base model: got: " << lp_base << ", expected: " << LP_BASE_STD);
	double dev = fobfilter.deviance_theta(THETA, PHI);
	BOOST_CHECK_MESSAGE((abs(dev - DEV_STD) < TEST_EPS),
		"Model deviance: got: " << dev << ", expected: " << DEV_STD);
}


BOOST_AUTO_TEST_CASE (hs_hd) {
	const char TSVNAME[] = "../sim-data/obrs_high-signal_high-damage_data.tsv";
	const double THETA_0 = 0.1;
	const double PHI_0 = 0.1;
	const double THETA_F_STD = .0784761;
	const double PVAL_F_STD = 1.4061715e-18;
	const double THETA_CA_STD = 0.0797567;
	const double PHI_CA_STD = 0.0952038;
	const double PVAL_CA_STD = 9.9433447e-19;
	BOOST_TEST_MESSAGE("\nRunning high-signal, high-damage test:");
	common_filter_math_test(TSVNAME, THETA_0, PHI_0, THETA_F_STD, PVAL_F_STD, THETA_CA_STD, PHI_CA_STD,
		PVAL_CA_STD);
}

BOOST_AUTO_TEST_CASE (hs_ld) {
	const char TSVNAME[] = "../sim-data/obrs_high-signal_low-damage_data.tsv";
	const double THETA_0 = 0.1;
	const double PHI_0 = 0.01;
	const double THETA_F_STD = .0879539;
	const double PVAL_F_STD = 2.6625633e-26;
	const double THETA_CA_STD = 0.1096559;
	const double PHI_CA_STD = 3.0590232e-7;
	const double PVAL_CA_STD = 2.1828736e-47;
	BOOST_TEST_MESSAGE("\nRunning high-signal, low-damage test:");
	common_filter_math_test(TSVNAME, THETA_0, PHI_0, THETA_F_STD, PVAL_F_STD, THETA_CA_STD, PHI_CA_STD,
		PVAL_CA_STD);
}

BOOST_AUTO_TEST_CASE ( ls_hd ) {
	const char TSVNAME[] = "../sim-data/obrs_low-signal_high-damage_data.tsv";
	const double THETA_0 = 0.01;
	const double PHI_0 = 0.1;
	const double THETA_F_STD = 3.0590232e-7;
	const double PVAL_F_STD = 1.0;
	const double THETA_CA_STD = 3.0590232e-7;
	const double PHI_CA_STD = .1223057;
	const double PVAL_CA_STD = 1.0;
	BOOST_TEST_MESSAGE("\nRunning low-signal, high-damage test:");
	common_filter_math_test(TSVNAME, THETA_0, PHI_0, THETA_F_STD, PVAL_F_STD, THETA_CA_STD, PHI_CA_STD,
		PVAL_CA_STD);
}

BOOST_AUTO_TEST_CASE (ls_ld) {
	const char TSVNAME[] = "../sim-data/obrs_low-signal_low-damage_data.tsv";
	const double THETA_0 = 0.01;
	const double PHI_0 = 0.01;
	const double THETA_F_STD = 3.0590231e-7;
	const double PVAL_F_STD = 1.0;
	const double THETA_CA_STD = 3.0590232e-7;
	const double PHI_CA_STD = .0088582;
	const double PVAL_CA_STD = 1.0;
	BOOST_TEST_MESSAGE("\nRunning low-signal, low-damage test:");
	common_filter_math_test(TSVNAME, THETA_0, PHI_0, THETA_F_STD, PVAL_F_STD, THETA_CA_STD, PHI_CA_STD,
		PVAL_CA_STD);
}

BOOST_AUTO_TEST_CASE (ms_hd) {
	const char TSVNAME[] = "../sim-data/obrs_med-signal_high-damage_data.tsv";
	const double THETA_0 = 0.05;
	const double PHI_0 = 0.1;
	const double THETA_F_STD = .0523259;
	const double PVAL_F_STD = 7.6402222e-9;
	const double THETA_CA_STD = .0416997;
	const double PHI_CA_STD = .141235;
	const double PVAL_CA_STD = 4.1888534e-8;
	// basic initialization
	BOOST_TEST_MESSAGE("\nRunning med-signal, high-damage test:");
	common_filter_math_test(TSVNAME, THETA_0, PHI_0, THETA_F_STD, PVAL_F_STD, THETA_CA_STD, PHI_CA_STD,
		PVAL_CA_STD);
}

BOOST_AUTO_TEST_CASE (ms_ld) {
	const char TSVNAME[] = "../sim-data/obrs_med-signal_low-damage_data.tsv";
	const double THETA_0 = 0.05;
	const double PHI_0 = 0.01;
	const double THETA_F_STD = .0249587;
	const double PVAL_F_STD = 1.4061715e-18;
	const double THETA_CA_STD = .0397191;
	const double PHI_CA_STD = 3.0590236e-7;
	const double PVAL_CA_STD = 9.9433447e-19;
	// basic initialization
	BOOST_TEST_MESSAGE("\nRunning med-signal, low-damage test:");
	common_filter_math_test(TSVNAME, THETA_0, PHI_0, THETA_F_STD, PVAL_F_STD, THETA_CA_STD, PHI_CA_STD,
		PVAL_CA_STD);
}

BOOST_AUTO_TEST_CASE (ns_hd) {
	const char TSVNAME[] = "../sim-data/obrs_no-signal_high-damage_data.tsv";
	const double THETA_0 = 0.0;
	const double PHI_0 = 0.1;
	const double THETA_F_STD = 3.0590232e-7;
	const double PVAL_F_STD = 1.0;
	const double THETA_CA_STD = 3.0590232e-7;
	const double PHI_CA_STD = .0620615;
	const double PVAL_CA_STD = 1.0;
	BOOST_TEST_MESSAGE("\nRunning no-signal, high-damage test:");
	common_filter_math_test(TSVNAME, THETA_0, PHI_0, THETA_F_STD, PVAL_F_STD, THETA_CA_STD, PHI_CA_STD,
		PVAL_CA_STD);
}

BOOST_AUTO_TEST_CASE (ns_ld) {
	const char TSVNAME[] = "../sim-data/obrs_no-signal_low-damage_data.tsv";
	const double THETA_0 = 0.0;
	const double PHI_0 = 0.01;
	const double THETA_F_STD = 3.0590232e-7;
	const double PVAL_F_STD = 1.0;
	const double THETA_CA_STD = 3.0590232e-7;
	const double PHI_CA_STD = .0093751;
	const double PVAL_CA_STD = 1.0;
	BOOST_TEST_MESSAGE("\nRunning no-signal, low-damage test:");
	common_filter_math_test(TSVNAME, THETA_0, PHI_0, THETA_F_STD, PVAL_F_STD, THETA_CA_STD, PHI_CA_STD,
		PVAL_CA_STD);
}

BOOST_AUTO_TEST_SUITE_END()
