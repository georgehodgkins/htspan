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
#include "htspan/io/faidx_reader.hpp"
#include "htspan/freq_orient_bias_quant.hpp"
#include "htspan/bayes_orient_bias_quant.hpp"

#include "htspan/io/simul_writer.hpp"
using namespace hts;

namespace hts {

/**
* This function provides an interface between the
* calling frontend and the backend classes for
* damage quantification.
*
* @param ref Reference nucleotide
* @param alternative nucleotide
* @param align_fname path to BAM file containing sequence of interest
* @param ref_fname path to FASTA file containing the corresponding reference sequence
* @param min_mapq Minimum mapping quality for a read to be considered
* @param min_baseq Minimum base quality for a read to be considered
* @param keep_dup Whether to exclude duplicated reads from consideration
* @param max_qreads Maximum number of (passing) reads to analyze
*/
bool orient_bias_quantify_freq(nuc_t ref, nuc_t alt, double min_mapq, double min_baseq, bool keep_dup,
		 size_t max_qreads, piler &p, faidx_reader &faidx, bool plain_output) {
	// Initialize objects
	freq_orient_bias_quant_f obquant (ref, alt);
	//set quality filter theresholds
	p.qfilter.min_mapq = min_mapq;
	p.qfilter.min_baseq = min_baseq;
	if (keep_dup) {
		p.qfilter.disable_excl_flags(BAM_FDUP);
	}
	
	//process reads up to max_reads
	size_t n_reads = 0;

	// iterator through pileup positions
	size_t j = 0;
	while (n_reads < max_qreads) {
		const bam_pileup1_t* pile = p.next();
		if (pile == NULL) break;

		size_t n = p.size();

		frontend::global_log.v(3) << "pileup " << j << ": " << n << " reads at tid = " 
			   << p.tid << ", pos = " << p.pos << '\n';

		const char* seq = faidx.get(p.tid, p.pos, p.pos+1);
		if (seq == NULL) {
			frontend::global_log.v(1) << "Warning: reference sequence could not be retrieved for contig " << p.tid
				<< " at position " << p.pos << '\n';
			continue;
		}
	
		
		// Check here whether reference nucleotide at position is
		//     G or C (oxoG: G>T, C>A)
		//     C or G (FFPE: C>T, G>A)
		char s_ref = char_to_nuc(seq[0]);
		frontend::global_log.v(3) << "ref = " << nuc_to_char(ref) << '\n';
		if (s_ref == ref || s_ref == nuc_complement(ref)) {
			// accumulate statistics
			n_reads += obquant.push(pile, n, p.pos);
		}
		++j;
	}
	//calculate the phi estimator
	double phi = obquant();
	//output the estimator
	if (!plain_output) {
		frontend::global_log.v(1) << "Phi estimator: ";
	} 
	frontend::global_log.v(1) << phi << '\n';
	return true;
}

}// namespace hts