#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <cstdlib>
using namespace std;


#include <htslib/hts.h>
#include <htslib/sam.h>

#include "fetcher.hpp"

int main(int argc, char** argv) {
	--argc;

	if (argc != 6) {
		cerr << "usage: " << argv[0]
			<< " <bam_path> <target> <pos> <keep_dup> <min_mapq> <min_baseq>" << endl;
		return 1;
	}

	char* path = argv[1];
	char* target = argv[2];
	int32_t pos = atol(argv[3]) - 1;
	bool keep_dup = (atoi(argv[4]) != 0);
	int min_mapq = atoi(argv[5]);
	int min_baseq = atoi(argv[6]);

	fetcher f;

	// disable quality filters
	f.qfilter.min_mapq = min_mapq;
	f.qfilter.min_baseq = min_baseq;

	if (keep_dup) {
		f.qfilter.disable_excl_flags(BAM_FDUP);
	}

	//char target[] = "chr17";
	//int32_t pos = 7674420 - 1;   // dup:   8 A, 94 G
	                               // nodup: 2 A, 24 G
	//int32_t pos = 7674360 - 1;   // dup:   0 A, 759 C, 0 G, 9 T
	                               // nodup: 0 A, 198 C, 0 G, 1 T
	//int32_t pos = 7674361 - 1;   // dup:   1 A, 2 C, 0 G, 536 T, 202 -, 15 +
	                               // nodup: 0 A, 0 C, 0 G, 154 T,  40 -,  4 +
	
	if (!f.open(path)) {
		cerr << "ERROR: could not open BAM file for reading: " << path << endl;
		return 1;
	}

	int32_t rid = bam_name2id(f.hdr, target);
	f.fetch(rid, pos);

	size_t n = f.pile.queries.size();
	size_t nA = 0, nC = 0, nG = 0, nT = 0, nN = 0, nNull = 0, nDel = 0, nIns = 0;
	for (size_t i = 0; i < n; ++i) {
		bam1_t *b = f.pile.queries[i];
		switch (query_nucleotide(b, pos)) {
		case nuc_A:
			++nA;
			break;
		case nuc_C:
			++nC;
			break;
		case nuc_G:
			++nG;
			break;
		case nuc_T:
			++nT;
			break;
		case nuc_N:
			++nN;
			break;
		case nuc_NULL:
			++nNull;
			break;
		case nuc_DEL:
			++nDel;
			break;
		case nuc_INS:
			++nIns;
			break;
		}
	}

	// print summaries
	cout << "A " << nA << endl;
	cout << "C " << nC << endl;
	cout << "G " << nG << endl;
	cout << "T " << nT << endl;
	cout << "N " << nN << endl;
	cout << "- " << nDel << endl;
	cout << "+ " << nIns << endl;
	cout << ". " << nNull << endl;
	cout << "total " << n << endl;

	return 0;
}
