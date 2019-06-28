#define BOOST_TEST_MODULE io_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <typeinfo>

#include "htspan/nucleotide.hpp"
#include "htspan/io/snv_reader.hpp"
#include "htspan/io/snv_writer.hpp"
#include "htspan/orient_bias_data.hpp"
#include "htspan/bayes_orient_bias_filter.hpp"

#include "test.hpp"

// NB: this test is designed to work with specific files: test-io.tsv/vcf in data dir
template <typename SnvReader>
void common_snvr_test (const char SNVNAME[]) {
	BOOST_TEST_MESSAGE("Starting SNV reader test for type " << typeid(SnvReader).name() << " and path " << SNVNAME);
	SnvReader snvr (SNVNAME);
	using hts::snv::record;
	const record L1 ("19", 110, 'A', 'C');
	const record L2A ("19", 111, 'A', 'G');
	const record L2B ("19", 111, 'A', 'C');
	const record L2C ("19", 111, 'A', 'T');
	const int L3_ERR = 1;
	const int L4_ERR = 1;
	const int L5_ERR = 3;
	const int L6_ERR = 3;
	const int L7_ERR = 3;
	record rec;
	// LINE 1: A>C
	snvr.next(rec);
	BOOST_CHECK_MESSAGE(rec == L1 && snvr.error() == 0,
		"Record does not match or read error was encountered." << 
		"\nGot: " << rec.to_string() <<
		"\nExp: " << L1.to_string() << 
		"\nErr: " << snvr.error());

	// LINE 2: A>G,C,T; should split alts
	snvr.next(rec);
	BOOST_CHECK_MESSAGE(rec == L2A && snvr.error() == 0,
		"Record does not match or read error was encountered." << 
		"\nGot: " << rec.to_string() <<
		"\nExp: " << L2A.to_string() << 
		"\nErr: " << snvr.error());
	snvr.next(rec);
	BOOST_CHECK_MESSAGE(rec == L2B && snvr.error() == 0,
		"Record does not match or read error was encountered." << 
		"\nGot: " << rec.to_string() <<
		"\nExp: " << L2B.to_string() << 
		"\nErr: " << snvr.error());
	snvr.next(rec);
	BOOST_CHECK_MESSAGE(rec == L2C && snvr.error() == 0,
		"Record does not match or read error was encountered." << 
		"\nGot: " << rec.to_string() <<
		"\nExp: " << L2C.to_string() << 
		"\nErr: " << snvr.error());

	// LINE 3: N>A; expect error 1 (zero-len allele)
	snvr.next(rec);
	BOOST_CHECK_MESSAGE(snvr.error() == L3_ERR,
		"Did not get expected error." <<
		"\nGot err: " << snvr.error() << " Exp err: " << L3_ERR <<
		"\nGot rec: " << rec.to_string());

	// LINE 4: T>N; expect error 1 (zero-len allele)
	snvr.next(rec);
	BOOST_CHECK_MESSAGE(snvr.error() == L4_ERR,
		"Did not get expected error." <<
		"\nGot err: " << snvr.error() << " Exp err: " << L4_ERR <<
		"\nGot rec: " << rec.to_string());
	// LINE 5: A>GT; expect error 3 (multi-len allele)
	snvr.next(rec);
	BOOST_CHECK_MESSAGE(snvr.error() == L5_ERR,
		"Did not get expected error." <<
		"\nGot err: " << snvr.error() << " Exp err: " << L5_ERR <<
		"\nGot rec: " << rec.to_string());
	// LINE 6: AT>G; expect error 3 (zero-len allele)
	snvr.next(rec);
	BOOST_CHECK_MESSAGE(snvr.error() == L6_ERR,
		"Did not get expected error." <<
		"\nGot err: " << snvr.error() << " Exp err: " << L6_ERR <<
		"\nGot rec: " << rec.to_string());
	// LINE 7: G>GA,GAC; expect error 3 (zero-len allele)
	snvr.next(rec);
	BOOST_CHECK_MESSAGE(snvr.error() == L7_ERR,
		"Did not get expected error." <<
		"\nGot err: " << snvr.error() << " Exp err: " << L7_ERR <<
		"\nGot rec: " << rec.to_string());
}

template <typename SnvReader, typename SnvWriter>
void common_snvw_test (const char IN_SNVNAME[], const char OUT_SNVNAME[]) {
	BOOST_TEST_MESSAGE("Starting SNV writer test for reader/writer types " << typeid(SnvReader).name() << '/' <<
		typeid(SnvWriter).name() << " and input path " << IN_SNVNAME);
	SnvReader snvr(IN_SNVNAME);
	SnvWriter snvw(OUT_SNVNAME, snvr);
	hts::orient_bias_data foo (nuc_G, nuc_T, 0);
	hts::bayes_orient_bias_filter_f bobfilter(foo);
	hts::snv::record rec, recF, recL;
	recL.pos = 1;
	snvw.add_filter_tag(bobfilter);
	snvr.next(recF);//discard first line
	snvr.next(recF); // store second line (new first line)
	do {
		if (recL.pos % 2) {
			snvw.write(snvr.get_underlying());
		} else {
			snvw.write(snvr.get_underlying(), bobfilter.text_id);
		}
	} while (snvr.next(recL)); // store last line
	snvw.close();
	snvr.close();
	snvr.open(OUT_SNVNAME);// read SNVs back in for checking
	snvr.next(rec);// get first record
	BOOST_CHECK_MESSAGE(rec.pos == recF.pos,
		"First record in written SNV does not match.");
	while (snvr.next(rec)) {}// get last record
	BOOST_CHECK_MESSAGE(rec.pos == recL.pos,
		"Last record in written SNV does not match.");
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

BOOST_AUTO_TEST_CASE (tsv_snvr) {
	const char SNVNAME[] = "../../data/test-io.tsv";
	common_snvr_test<hts::snv::tsv_reader>(SNVNAME);

}

BOOST_AUTO_TEST_CASE (uncompressed_vcf_snvr) {
	const char SNVNAME[] = "../../data/test-io.vcf";
	common_snvr_test<hts::snv::vcf_reader>(SNVNAME);
}

BOOST_AUTO_TEST_CASE (gzip_vcf_snvr) {
	const char SNVNAME[] = "../../data/test-io.vcf.gz";
	common_snvr_test<hts::snv::vcf_reader>(SNVNAME);
}

BOOST_AUTO_TEST_CASE (bgzip_vcf_snvr) {
	const char SNVNAME[] = "../../data/test-io.vcf.bgz";
	common_snvr_test<hts::snv::vcf_reader>(SNVNAME);
}

BOOST_AUTO_TEST_CASE (tsv_snvw) {
	const char IN_SNVNAME[] = "../../data/snv.tsv";
	const char OUT_SNVNAME[] = "../../data/test_out.tsv";
	common_snvw_test<hts::snv::tsv_reader, hts::snv::tsv_writer>(IN_SNVNAME, OUT_SNVNAME);
}

BOOST_AUTO_TEST_CASE (uncompressed_vcf_snvw) {
	const char IN_SNVNAME[] = "../../data/sample.vcf";
	const char OUT_SNVNAME[] = "../../data/test_out.vcf";
	common_snvw_test<hts::snv::vcf_reader, hts::snv::vcf_writer>(IN_SNVNAME, OUT_SNVNAME);
}

BOOST_AUTO_TEST_CASE (gzip_vcf_snvw) {
	const char IN_SNVNAME[] = "../../data/sample.vcf.gz";
	const char OUT_SNVNAME[] = "../../data/test_out.vcf.gz";
	common_snvw_test<hts::snv::vcf_reader, hts::snv::vcf_writer>(IN_SNVNAME, OUT_SNVNAME);
}

BOOST_AUTO_TEST_CASE (bgzip_vcf_snvw) {
	const char IN_SNVNAME[] = "../../data/sample.vcf.bgz";
	const char OUT_SNVNAME[] = "../../data/test_out.vcf.bgz";
	common_snvw_test<hts::snv::vcf_reader, hts::snv::vcf_writer>(IN_SNVNAME, OUT_SNVNAME);
}

BOOST_AUTO_TEST_SUITE_END()