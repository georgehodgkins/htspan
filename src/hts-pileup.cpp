#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <cstdlib>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "htspan/piler.hpp"
#include "htspan/bam.hpp"

using namespace hts;

int main(int argc, char** argv) {
	--argc;

	if (argc < 4) {
		cerr << "usage: " << argv[0]
			<< " <bam_path> <keep_dup> <min_mapq> <min_baseq> [<n_pileups:5>]" << endl;
		return 1;
	}

	char* path = argv[1];
	bool keep_dup = (atoi(argv[2]) != 0);
	int min_mapq = atoi(argv[3]);
	int min_baseq = atoi(argv[4]);
	size_t n_pileups = 5;
	if (argc > 4) {
		n_pileups = atoi(argv[5]);
	}

	piler p;
	
	// set quality filters
	p.qfilter.min_mapq = min_mapq;
	p.qfilter.min_baseq = min_baseq;

	if (keep_dup) {
		p.qfilter.disable_excl_flags(BAM_FDUP);
	}

	if (!p.open(path)) {
		cerr << "ERROR: could not open BAM file for reading: " << path << endl;
		return 1;
	}

	// output header
	cout << "read#\tqpos\tnuc\tbaseq\tmapq\tqname\n";

	// iterate through pileup positions, ignoring those with no passing reads
	for (size_t j = 0; j < n_pileups; ++j) {
		const vector<bam1_t*> &pile = p.next();
		if (pile.empty()) {
			--j;
			continue;
		}
		size_t n = p.size();

		cout << "pileup " << j << ": " << n << " passing reads at tid = " << p.tid << ", pos = " << p.pos << endl;

		for (size_t i = 0; i < n; ++i) {
			const uint8_t* seq = bam_get_seq(pile[i]);
			const uint8_t* qual = bam_get_qual(pile[i]);
			int32_t qpos = query_position(pile[i], p.pos);
			cout
				<< i << '\t'
				<< qpos << '\t'
				<< nuc_to_char(bam_seqi(seq, qpos)) << '\t'
				<< (int) qual[qpos] << '\t' 
				<< (int) pile[i]->core.qual << '\t'
				<< bam_get_qname(pile[i]) << endl;
		}
	}

	return 0;
}
