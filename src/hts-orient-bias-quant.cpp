#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <cstdlib>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "htspan/piler.hpp"
#include "htspan/faidx_reader.hpp"
#include "htspan/orient_bias_quant.hpp"
using namespace hts;

int main(int argc, char** argv) {
	--argc;

	if (argc != 6) {
		cerr << "usage: " << argv[0]
			<< " <bam_path> <ref_fasta> <max_reads> <keep_dup> <min_mapq> <min_baseq>" << endl;
		return 1;
	}

	char* path = argv[1];
	char* ref_fasta = argv[2];
	int max_reads = atoi(argv[3]);
	bool keep_dup = (atoi(argv[4]) != 0);
	int min_mapq = atoi(argv[5]);
	int min_baseq = atoi(argv[6]);

	piler p;
	faidx_reader faidx;
	
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

	// open reference genome
	if (!faidx.open(ref_fasta)) {
		cerr << "ERROR: could not open FASTA file for reading: " << path << endl;
		return 1;
	}
	
	// oxoG: G>T on read 1
	// FFPE: C>T on read 1
	orient_bias_quant_f obquant(nuc_G, nuc_T);

	size_t n_reads = 0;

	// iterator through pileup positions
	size_t j = 0;
	while (n_reads < max_reads) {
		const bam_pileup1_t* pile = p.next();
		if (pile == NULL) break;

		size_t n = p.size();

		cerr << "pileup " << j << ": " << n << " reads at tid = " 
			   << p.tid << ", pos = " << p.pos << endl;

		const char* seq = faidx.get(p.tid, p.pos, p.pos+1);
		if (seq == NULL) {
			cerr << "WARNING: reference sequence could not be retrieved for contig " << p.tid
				<< " at position " << p.pos << endl;
			continue;
		}
	
		
		// Check here whether reference nucleotide at position is
		//     G or C (oxoG: G>T, C>A)
		//     C or G (FFPE: C>T, G>A)
		char ref = char_to_nuc(seq[0]);
		cerr << "ref = " << nuc_to_char(ref) << endl;
		if (ref == nuc_G || ref == nuc_C) {
			// accumulate statistics
			n_reads += obquant.push(pile, n, p.pos);
		}

		++j;
	}

	// estimate of global damage
	double phi = obquant();

	cerr << "xi = " << obquant.xi << ", ni = " << obquant.ni << endl;
	cerr << "xc = " << obquant.xc << ", nc = " << obquant.nc << endl;
	
	cout << "phi = " << phi << endl;

	return 0;
}
