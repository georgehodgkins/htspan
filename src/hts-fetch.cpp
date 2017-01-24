#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "fetcher.hpp"

int main(int argc, char** argv) {
	char path[] = "../data/KP7-092916.bam"; 

	fetcher f;

	f.open(path);

	// set target region
	char target[] = "chr17";
	int32_t rid = bam_name2id(f.hdr, target);
	int32_t pos = 7674420 - 1;  // 8 A, 94 G
	//int32_t pos = 7674361 - 1;   // 1 A, 2 C, 0 G, 536 T, 202 del, 15 ins
	//int32_t pos = 7674360 - 1;    // 0 A, 759 C, 0 G, 9 T
	
	//f.fetch(rid, pos);
	f.fetch(rid, pos, nuc_A, true);
	//f.fetch(rid, pos, nuc_C, true);

	for (size_t i = 0; i < f.pile.queries.size(); ++i) {
		bam1_t *b = f.pile.queries[i];
		//print_query_fastq(b);
		print_query_fasta(b);
	}	
	cout << endl;
	for (size_t i = 0; i < f.pile.mates.size(); ++i) {
		bam1_t *b = f.pile.mates[i];
		//print_query_fastq(b);
		print_query_fasta(b);
	}	

	return 0;
}
