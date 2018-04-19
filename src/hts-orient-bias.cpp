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


int main(int argc, char** argv) {
	--argc;

	if (argc != 7) {
		cerr << "usage: " << argv[0]
			<< " <bam_path> <target> <pos> <keep_dup> <min_mapq> <min_baseq> <oxidation_error>" << endl;
		return 1;
	}

	char* path = argv[1];
	char* target = argv[2];
	int32_t pos = atol(argv[3]) - 1;
	bool keep_dup = (atoi(argv[4]) != 0);
	int min_mapq = atoi(argv[5]);
	int min_baseq = atoi(argv[6]);
	double phi = atof(argv[7]);

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

	size_t n = f.pile.queries.size();

	orient_bias_filter_f obfilter(nuc_G, nuc_T, n);

	size_t nA = 0, nC = 0, nG = 0, nT = 0, nN = 0, nNull = 0, nDel = 0, nIns = 0;
	for (size_t i = 0; i < n; ++i) {
		bam1_t *b = f.pile.queries[i];
		obfilter.push(b, pos);
	}  // for (size_t i = 0; i < n; ++i)

	cout << "# theta_hat = " << obfilter.estimate_theta_given(phi) << endl;
	cout << "# p = " << obfilter(phi) << endl;

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

	return 0;
}
