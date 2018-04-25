#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <cstdlib>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "htspan/piler.hpp"
using namespace hts;

int main(int argc, char** argv) {
	--argc;

	if (argc != 4) {
		cerr << "usage: " << argv[0]
			<< " <bam_path> <keep_dup> <min_mapq> <min_baseq>" << endl;
		return 1;
	}

	char* path = argv[1];
	bool keep_dup = (atoi(argv[2]) != 0);
	int min_mapq = atoi(argv[3]);
	int min_baseq = atoi(argv[4]);

	piler p;
	
	// disable quality filters
	p.qfilter.min_mapq = min_mapq;
	p.qfilter.min_baseq = min_baseq;

	if (keep_dup) {
		p.qfilter.disable_excl_flags(BAM_FDUP);
	}

	if (!p.open(path)) {
		cerr << "ERROR: could not open BAM file for reading: " << path << endl;
		return 1;
	}

	// iterator through pileup positions
	for (size_t j = 0; j < 5; ++j) {
		const bam_pileup1_t* pile = p.next();
		size_t n = p.size();

		cout << "pileup " << j << ": " << n << " reads at tid = " << p.tid << ", pos = " << p.pos << endl;

		for (size_t i = 0; i < n; ++i) {
			int32_t qpos = pile[i].qpos;
			const bam1_t* b = pile[i].b;
			const uint8_t* seq = bam_get_seq(b);
			const uint8_t* qual = bam_get_qual(b);
			cout
				<< i << '\t'
				<< bam_get_qname(b) << '\t'
				<< qpos << '\t'
				<< nuc_to_char(bam_seqi(seq, qpos)) << '\t'
				<< (int) qual[qpos] << '\t' 
				<< (int) b->core.qual << endl;
		}
	}

	return 0;
}
