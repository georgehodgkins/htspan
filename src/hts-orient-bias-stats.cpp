#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <cstdlib>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "htspan/fetcher.hpp"
#include "htspan/orient_bias_filter.hpp"
using namespace hts;

void print_orient_bias_filter(const orient_bias_filter_f& obfilter) {
	cout << "base\torient\terror\tcnuc\tmember\tstrand\tqual" << endl;
	for (size_t i = 0; i < obfilter.size(); ++i) {
		cout
			<< (int) obfilter.bases[i] << '\t'
			<< obfilter.orients[i] << '\t'
			<< obfilter.errors[i] << '\t'
			<< obfilter.cnucs[i] << '\t'
			<< obfilter.members[i] << '\t'
			<< obfilter.strands[i] << '\t'
			<< obfilter.quals[i] << endl;
	}
}

int main(int argc, char** argv) {
	--argc;

	if (argc != 9) {
		cerr << "usage: " << argv[0]
			<< " <bam_path> <target> <pos> <ref> <alt> <phi> <keep_dup> <min_mapq> <min_baseq>" << endl;
		return 1;
	}

	char* path = argv[1];
	char* target = argv[2];
	int32_t pos = atol(argv[3]) - 1;
	nuc_t nt_ref = char_to_nuc(argv[4][0]);
	nuc_t nt_alt = char_to_nuc(argv[5][0]);
	double phi = atof(argv[6]);
	bool keep_dup = (atoi(argv[7]) != 0);
	int min_mapq = atoi(argv[8]);
	int min_baseq = atoi(argv[9]);

	gsl_set_error_handler_off();

	fetcher f;

	// disable quality filters
	f.qfilter.min_mapq = min_mapq;
	f.qfilter.min_baseq = min_baseq;

	if (keep_dup) {
		f.qfilter.disable_excl_flags(BAM_FDUP);
	}
	
	if (!f.open(path)) {
		cerr << "ERROR: could not open BAM file for reading: " << path << endl;
		return 1;
	}

	// fetch reads at target position
	int32_t rid = bam_name2id(f.hdr, target);
	f.fetch(rid, pos);

	// test for oxoG damage (G>T on read 1 or C>A on read 2)
	orient_bias_filter_f obfilter(nuc_G, nuc_T, f.pile.queries.size());
	obfilter.push(f.pile.queries, pos, nt_ref, nt_alt);

	cout << "# Variant test adjusted for oxoG damage" << endl;
	cout << "# theta_hat = " << obfilter.estimate_theta_given(phi, obfilter.estimate_initial_theta()) << endl;
	cout << "# p = " << obfilter(phi) << endl;

	print_orient_bias_filter(obfilter);

	return 0;
}
