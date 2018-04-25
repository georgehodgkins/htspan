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

	
	query_filter_f qfilter;

	// disable quality filters
	qfilter.min_mapq = min_mapq;
	qfilter.min_baseq = min_baseq;

	if (keep_dup) {
		qfilter.disable_excl_flags(BAM_FDUP);
	}

	//char target[] = "chr17";
	//int32_t pos = 7674420 - 1;   // dup:   8 A, 94 G
	                               // nodup: 2 A, 24 G
	//int32_t pos = 7674360 - 1;   // dup:   0 A, 759 C, 0 G, 9 T
	                               // nodup: 0 A, 198 C, 0 G, 1 T
	//int32_t pos = 7674361 - 1;   // dup:   1 A, 2 C, 0 G, 536 T, 202 -, 15 +
	                               // nodup: 0 A, 0 C, 0 G, 154 T,  40 -,  4 +
	
	piler p;

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
			uint8_t* seq = bam_get_seq(pile[i].b);
			cout << i << '\t' << bam_get_qname(pile[i].b) << '\t' << 
				qpos << '\t' << nuc_to_char(bam_seqi(seq, qpos)) << endl;
		}
	}

	return 0;
}
