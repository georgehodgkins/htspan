#include <iostream>
#include <cstring>
#include <cstdlib>
using namespace std;

#include "htspan/fetcher.hpp"
using namespace hts;

int main(int argc, char** argv) {
	--argc;

	if (argc != 3 && argc != 7) {
		cerr << "usage: " << argv[0] 
			<< "<bam_path> <target> <pos> [<nucleotide> <keep_dup> <get_mate> <format>]" << endl;
		return 1;
	}

	char* path = argv[1];
	char* target = argv[2];
	int32_t pos = atol(argv[3]) - 1;

	uint8_t nucleotide = nuc_N;
	bool keep_dup = false;
	bool get_mate = false;
	bool use_fastq = false;

	if (argc == 7) {
		nucleotide = char_to_nuc(argv[4][0]);
		keep_dup = (atoi(argv[5]) != 0);
		get_mate = (atoi(argv[6]) != 0);
		char* format = argv[7];
		if (std::strcmp(format, "fastq") == 0 || std::strcmp(format, "fq") == 0) {
			use_fastq = true;
		}
	}

	cerr << path << ':' << target << ':' << pos << ':' << nuc_to_char(nucleotide) << endl;
	
	fetcher f;

	if (keep_dup) {
		f.qfilter.disable_excl_flags(BAM_FDUP);
	}

	if (!f.open(path)) {
		cerr << "ERROR: could not open BAM file for reading: " << path << endl;
		return 1;
	}
	
	int32_t rid = bam_name2id(f.hdr, target);
	f.fetch(rid, pos, nucleotide, get_mate);

	if (get_mate) {
		// get mate reads as well
		if (use_fastq) {
			for (size_t i = 0; i < f.pile.queries.size(); ++i) {
				bam1_t *b1 = f.pile.queries[i];
				print_query_fastq(b1);
				bam1_t *b2 = f.pile.mates[i];
				print_query_fastq(b2);
				cout << endl;
			}	
		} else {
			for (size_t i = 0; i < f.pile.queries.size(); ++i) {
				bam1_t *b1 = f.pile.queries[i];
				print_query_fasta(b1);
				bam1_t *b2 = f.pile.mates[i];
				print_query_fasta(b2);
				cout << endl;
			}
		}
	} else {
		// only get reads mapping to position
		if (use_fastq) {
			for (size_t i = 0; i < f.pile.queries.size(); ++i) {
				bam1_t *b1 = f.pile.queries[i];
				print_query_fastq(b1);
			}	
		} else {
			for (size_t i = 0; i < f.pile.queries.size(); ++i) {
				bam1_t *b1 = f.pile.queries[i];
				print_query_fasta(b1);
			}
		}
	}

	return 0;
}
