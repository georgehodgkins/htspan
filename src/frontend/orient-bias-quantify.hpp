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

#include "simul_writer.hpp"
using namespace hts;

/**
* This file contains driver functions that implement the quantification process.
*
* Alignment data from an hts::piler is compared against an indexed reference in
* an hts::faidx by an hts::<model>_orient_bias_quant_f object, which performs analysis based
* on the set of observed nucleotides at each site compared to the reference
* and returns estimates of parameters used in the identification process.
*/

namespace hts {

/**
* This function provides a driver for the process of frequentist quantification.
*
* Verbosity output levels:
*	1: Warnings, Phi estimator
*	2: Nothing extra
*	3: Info dump for each piled locus
*
* @param ref Reference nucleotide for damage type
* @param alt Alternative nucleotide for damage type
* @param min_mapq Minimum mapping quality for a read to be considered
* @param min_baseq Minimum base quality for a read to be considered
* @param keep_dup Whether to exclude duplicated reads from consideration
* @param max_qreads Maximum number of reads (that pass initial filters) to analyze
* @param p Already initialized hts::piler containing alignment (BAM data) of interest
* @param faidx Already intialized hts::faidx_reader containing reference sequence to compare against
* @param plain_output Whether to include non machine-readable output
* @return Whether the operation succeeded
*/
bool orient_bias_quantify_freq(nuc_t ref, nuc_t alt, double min_mapq, double min_baseq, bool keep_dup,
		 size_t max_qreads, piler &p, faidx_reader &faidx, bool plain_output) {
	// Initialize quantification class
	freq_orient_bias_quant_f fobquant (ref, alt);
	
	// set quality filter thresholds
	p.qfilter.min_mapq = min_mapq;
	p.qfilter.min_baseq = min_baseq;
	if (keep_dup) {
		p.qfilter.disable_excl_flags(BAM_FDUP);
	}
	
	// number of processed reads
	size_t n_reads = 0;

	// iterate through pileup positions until enough reads have been processed
	size_t j = 0;
	while (n_reads < max_qreads) {
		
		// advance to next locus
		const bam_pileup1_t* pile = p.next();
		if (pile == NULL) break;
		
		// number of reads at the locus
		size_t n = p.size();

		frontend::global_log.v(3) << "pileup " << j << ": " << n << " reads at tid = " 
			   << p.tid << ", pos = " << p.pos << '\n';

		// Get the reference sequence at the corresponding location
		const char* seq = faidx.get(p.tid, p.pos, p.pos+1);
		if (seq == NULL) {
			frontend::global_log.v(1) << "Warning: reference sequence could not be retrieved for contig " << p.tid
				<< " at position " << p.pos << '\n';
			continue;
		}
		
		// Check here whether reference nucleotide at position is damage-consistent
		char s_ref = char_to_nuc(seq[0]);
		frontend::global_log.v(3) << "ref = " << nuc_to_char(ref) << '\n';
		if (s_ref == ref || s_ref == nuc_complement(ref)) {
			// accumulate statistics
			n_reads += fobquant.push(pile, n, p.pos);
		}
		++j;
	}
	
	//calculate the phi estimator
	double phi = fobquant();
	
	//output the estimator
	if (!plain_output) {
		frontend::global_log.v(1) << "Phi estimator: ";
	} 
	frontend::global_log.v(1) << phi << '\n';
	return true;
} // end frequentist quant

// Compile-time fixed parameters for the stochastic gradient (stograd)
// optimization process used in Bayesian quantification

// Programatically, these do not need to be compile-time fixed,
// but it is not necessary to expose them to the end user
// Feel free to tweak them to achieve a better estimate

#define BATCH_SIZE 250 // size of each stograd batch
#define N_EPOCHS 5000 // maximum number of stograd epochs (stops early if convergence is achieved)
#define LRATE 1e-4 // stograd learning rate (scales size of each parameter adjustment)
#define CONVERGENCE_EPS 1e-6 // threshold for convergence (TODO: use command-line eps)

/**
* This function provides a driver for the process of Bayesian quantification.
*
* Verbosity output levels:
*	1: Warnings, alpha_phi and beta_phi estimates
*	2: alpha_theta and beta_theta estimates, E[phi] and E[theta]
*	3: Info dump for each piled locus
*
* @param ref Reference nucleotide for damage type
* @param alt Alternative nucleotide for damage type
* @param min_mapq Minimum mapping quality for a read to be considered
* @param min_baseq Minimum base quality for a read to be considered
* @param keep_dup Whether to exclude duplicated reads from consideration
* @param max_qreads Maximum number of reads (that pass initial filters) to analyze
* @param p Already initialized hts::piler containing alignment (BAM data) of interest
* @param faidx Already intialized hts::faidx_reader containing reference sequence to compare against
* @param plain_output Whether to include non machine-readable output
* @return Whether the operation succeeded
*/
bool orient_bias_quantify_bayes(nuc_t ref, nuc_t alt, double min_mapq, double min_baseq, bool keep_dup,
		 size_t max_qreads, piler &p, faidx_reader &faidx, bool plain_output) {
	// Initialize objects
	bayes_orient_bias_quant_f bobquant (ref, alt);
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
		
		// advance to next locus
		const bam_pileup1_t* pile = p.next();
		if (pile == NULL) break;
		
		// number of reads at the locus
		size_t n = p.size();

		frontend::global_log.v(3) << "pileup " << j << ": " << n << " reads at tid = " 
			   << p.tid << ", pos = " << p.pos << '\n';

		// Get reference sequence at that locus
		const char* seq = faidx.get(p.tid, p.pos, p.pos+1);
		if (seq == NULL) {
			frontend::global_log.v(1) << "Warning: reference sequence could not be retrieved for contig " << p.tid
				<< " at position " << p.pos << '\n';
			continue;
		}
	
		
		// Check whether reference nucleotide at position is damage-consistent
		char s_ref = char_to_nuc(seq[0]);
		frontend::global_log.v(3) << "ref = " << nuc_to_char(ref) << '\n';
		if (s_ref == ref || s_ref == nuc_complement(ref)) {
			// accumulate statistics
			n_reads += bobquant.push(pile, n, p.pos);
		}
		++j;
	}
	
	// Get hyperparameter estimates (computationally expensive)
	hts::hparams result = bobquant.operator()<stograd::stepper::constant<double> > (BATCH_SIZE, N_EPOCHS, LRATE, CONVERGENCE_EPS);

	// Output theta hyperparameters and E[theta], if verbosity >= 2
	if (!plain_output) {
		frontend::global_log.v(2) << "Alpha_theta: ";
	}
	frontend::global_log.v(2) << result.alpha_theta << '\t';
	if (!plain_output) {
		frontend::global_log.v(2) << "Beta_theta: ";
	}
	frontend::global_log.v(2) << result.beta_theta << '\t';
	if (!plain_output) {
		frontend::global_log.v(2) << "E[theta]: ";
	}
	frontend::global_log.v(2) << result.alpha_theta/(result.alpha_theta + result.beta_theta) << '\n';

	// Output phi hyperparameters, and E[phi] if verbosity >= 2
	if (!plain_output) {
		frontend::global_log.v(1) << "Alpha_phi: ";
	}
	frontend::global_log.v(1) << result.alpha_phi << '\t';
	if (!plain_output) {
		frontend::global_log.v(1) << "Beta_phi: ";
	}
	frontend::global_log.v(1) << result.beta_phi;
	frontend::global_log.v(2) << '\t';
	if (!plain_output) {
		frontend::global_log.v(2) << "E[phi]: ";
	}
	frontend::global_log.v(2) << result.alpha_phi/(result.alpha_phi + result.beta_phi);
	frontend::global_log.v(1) << '\n';
	
	return true;
} // end Bayesian quant

}// namespace hts