#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <cstdlib>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "htspan/fetcher.hpp"
#include "htspan/freq_orient_bias_filter.hpp"
#include "htspan/orient_bias_data.hpp"
using namespace hts;

void print_orient_bias_data(const orient_bias_data& data) {
	cout << "base\torient\terror\tcnuc\tmember\tstrand\tqual" << endl;
	for (size_t i = 0; i < data.size(); ++i) {
		cout
			<< (int) data.bases[i] << '\t'
			<< data.orients[i] << '\t'
			<< data.errors[i] << '\t'
			<< data.cnucs[i] << '\t'
			<< data.members[i] << '\t'
			<< data.strands[i] << '\t'
			<< data.quals[i] << endl;
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

	// process fetched data
	orient_bias_data data(nuc_G, nuc_T, f.pile.queries.size());
	data.push(f.pile.queries, pos, nt_ref, nt_alt);

	// test for oxoG damage (G>T on read 1 or C>A on read 2)
	freq_orient_bias_filter_f fobfilter(data);
	cout << "# Variant test adjusted for oxoG damage" << endl;
	cout << "# theta_hat = " << fobfilter.estimate_theta_given(phi) << endl;
	cout << "# p = " << fobfilter(phi) << endl;

	print_orient_bias_data(data);

	return 0;
}
