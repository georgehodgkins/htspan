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
#include "htspan/io/snv.hpp"
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

	if (argc != 6) {
		cerr << "usage: " << argv[0]
			<< " <bam_path> <snv.tsv> <phi> <keep_dup> <min_mapq> <min_baseq>" << endl;
		return 1;
	}

	const char* bam_path = argv[1];
	const char* snv_path = argv[2];
	double phi = atof(argv[3]);
	bool keep_dup = (atoi(argv[4]) != 0);
	int min_mapq = atoi(argv[5]);
	int min_baseq = atoi(argv[6]);

	gsl_set_error_handler_off();

	snv::reader snvr(snv_path);

	fetcher f;

	f.qfilter.min_mapq = min_mapq;
	f.qfilter.min_baseq = min_baseq;

	if (keep_dup) {
		f.qfilter.disable_excl_flags(BAM_FDUP);
	}
	
	if (!f.open(bam_path)) {
		cerr << "Error: could not open BAM file for reading: " << bam_path << endl;
		return 1;
	}

	cout << "p\ttheta_hat" << endl;
	snv::record rec;
	while (snvr.next(rec)) {
		// test for oxoG damage (G>T on read 1 or C>A on read 2)
		if ((rec.nt_ref == nuc_G && rec.nt_alt == nuc_T) ||
				(rec.nt_ref == nuc_C && rec.nt_alt == nuc_A)) {

			// fetch reads at target position
			int32_t rid = bam_name_to_id(f.hdr, rec.chrom);
			if (rid == -1) {
				cerr << "Warning: could not find " << rec.chrom << endl;
				continue;
			}

			if (!f.fetch(rid, rec.pos)) {
				cerr << "Warning: could not fetch reads for: " << rec.chrom << ':' << rec.pos << endl;
				continue;
			}

			cerr << "Info: fetched " << f.pile.queries.size() << " reads" << endl;

			// run filter
			orient_bias_filter_f obfilter(nuc_G, nuc_T, f.pile.queries.size());
			obfilter.push(f.pile.queries, rec.pos, rec.nt_ref, rec.nt_alt);

			cout << obfilter(phi) << '\t'
				<< obfilter.estimate_theta_given(phi) << endl;

			f.clear();
		} else {
			cout << "NA\tNA" << endl;
		}
	}

	return 0;
}
