#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <string>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "fetcher.hpp"
#include "mlat.hpp"

void print_align(const gfAlign& align) {
		cout << align.qStart << '\t' << align.qEnd << endl;
		cout << align.tStart << '\t' << align.tEnd << endl;
		cout << align.score << '\t' << align.matchCount << '\t' << align.mismatchCount << '\t' << align.repMatchCount << '\t' << align.qInsertCount << '\t' << align.qInsertBaseCount << '\t' << align.tInsertCount << '\t' << align.tInsertBaseCount << endl;
}

int main(int argc, char** argv) {
	char path[] = "../data/test.bam"; 

	fetcher f;

	f.open(path);

	// set target region
	char target[] = "chr17";
	int32_t rid = bam_name2id(f.hdr, target);
	int32_t pos = 7674420 - 1;  // dup: 8 A, 94 G; nodup: 2 A
	//int32_t pos = 7674361 - 1;   // 1 A, 2 C, 0 G, 536 T, 202 del, 15 ins
	//int32_t pos = 7674360 - 1;    // 0 A, 759 C, 0 G, 9 T
	
	//f.fetch(rid, pos);
	f.fetch(rid, pos, nuc_A, true);
	//f.fetch(rid, pos, nuc_C, true);
	

	mlat::Database m = mlat::Database("../data/TP53_hg38.2bit");

	for (size_t i = 0; i < f.pile.queries.size(); ++i) {
		bam1_t *b1 = f.pile.queries[i];
		bam1_t *b2 = f.pile.mates[i];

		string b1seq, b2seq;
		bam_seq_str(b1, b1seq);
		bam_seq_str(b2, b2seq);

		char querySeq[100];

		cout << i << endl;
		cout << bam_get_qname(b1) << endl;
		cout << b1seq << endl;
		strcpy(querySeq, b1seq.c_str());
		mlat::Result res1(m.search(querySeq));
		print_align(res1[0]);
		cout << res1[0].score - res1[1].score << endl;
		cout << bam_get_qname(b2) << endl;
		cout << b2seq << endl;
		strcpy(querySeq, b2seq.c_str());
		mlat::Result res2(m.search(querySeq));
		print_align(res2[0]);
		cout << res2[0].score - res2[1].score << endl;
		cout << res2[0].tEnd - res1[0].tStart << endl;

		cout << endl;


	}	

	return 0;
}
