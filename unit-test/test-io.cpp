#define BOOST_TEST_MODULE io_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

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

BOOST_AUTO_TEST_CASE (vcf_reader) {// UNDER CONSTRUCTION
	const char SNVNAME[] = "../../data/sample.vcf";
	const int RID_F = 20;
	const long int POS_F = 14370;
	const nuc_t REF_F = nuc_G;
	const nuc_t ALT_F = nuc_A;
	const int RID_L = 20;
	const long int POS_L = 1234567;
	const nuc_T REF_L = C;
	const nuc_T ALT_L = G;


BOOST_AUTO_TEST_SUITE_END()