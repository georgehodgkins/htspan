#include <iostream>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <numeric>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "htspan/piler.hpp"
#include "htspan/io/faidx_reader.hpp"
#include "htspan/freq_orient_bias_quant.hpp"
#include "htspan/bayes_orient_bias_quant.hpp"

#include "simul_writer.hpp"
#include "json.h"
using namespace hts;

/**
* This file contains driver functions that implement the quantification process.
*
* Alignment data from an hts::piler is compared against an indexed reference in
* an hts::faidx by an hts::<model>_orient_bias_quant_f object, which performs analysis based
* on the set of observed nucleotides at each site compared to the reference
* and returns estimates of parameters used in the identification process.
*
* This file is intentionally not include guarded; it contains driver code only and should not be
* included from multiple locations.
*/

// TODO: Redoc

namespace hts {


size_t obquant_accumulate (base_orient_bias_quant_f &obquant, piler &p, faidx_reader &faidx,
		size_t max_reads, bool plain_output) {

	// number of processed reads
	size_t n_reads = 0;

	// number of total reads (incl failing and irrelevant)
	size_t t_reads = 0;

	// number of processed sites
	size_t j = 0;

	// table header for info dump (only at verbosity==3)
	frontend::global_log.v(3) << "pileup\treads\ttid\tpos\tref\talts\txij\txcj\tnij\tncj\n";
	while (n_reads < max_reads) {
		
		// advance to next locus, break if EOF has been reached
		const vector<bam1_t*> &pile = p.next();
		if (p.n <= 0) break;
		
		// total number of reads at the locus, including those that failed the filter
		size_t n = p.size();
		t_reads += n;

		// Get the reference sequence at the corresponding location
		const char* seq = faidx.get(p.tid, p.pos, p.pos+1);
		if (seq == NULL) {
			frontend::global_log.v(1) << "Warning: reference sequence could not be retrieved for contig " << p.tid
				<< " at position " << p.pos << '\n';
			continue;
		}

		frontend::global_log.v(3) << j << '\t' << pile.size() << " (" << n << ")\t" << p.tid << '\t' << p.pos <<
			" (r" << (bam_is_read1(pile[0]) ? 1 : 2) << (bam_is_rev(pile[0]) ? '-' : '+') <<  ")\t" << seq[0] << '\t';

		// dump locus info and observed variables at verbosity==3
		// this adds some extra overhead to evaluate, plus the overhead inherent from large amounts of tty output
		if (frontend::global_log.v() >= 3 && !pile.empty()) {
			char palt;
			for (size_t i = 0; i < pile.size()-1; ++i) {
				palt = nuc_to_char(query_nucleotide(pile[i], p.pos));
				if (palt != seq[0]) {
					frontend::global_log.v(3) << palt << ',';
				}
			}
			palt = nuc_to_char(query_nucleotide(pile[pile.size()-1], p.pos));
			if (palt != seq[0]) {
				frontend::global_log.v(3) << palt;
			}
			frontend::global_log.v(3) << '\t';
		}

		// Check here whether reference nucleotide at position is damage-consistent
		nuc_t s_ref = char_to_nuc(seq[0]);
		
		// accumulate statistics
		n_reads += obquant.push(pile, p.pos, s_ref);
		//frontend::global_log.v(3) << obquant.xij() << '\t' << obquant.xcj() << '\t'
		//	<< obquant.nij() << '\t' << obquant.ncj() << '\t' << n_reads << " (" << t_reads << ")";
		frontend::global_log.v(3) << '\n';
		++j;
	}

	frontend::global_log.v(2) << "Reads processed: " << t_reads << " (passing: " << n_reads << "), sites: " << j << '\n';

	return n_reads;
}

/**
* This function provides a driver for the process of frequentist quantification.
*
* Verbosity output levels:
*	1: Warnings, Phi estimator
*	2: Overall observed variables, theta estimator
*	3: Info dump for each piled locus (pos, read count, observed vars)
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
bool orient_bias_quantify_freq(nuc_t ref, nuc_t alt, piler &p, faidx_reader &faidx, simpleson::json::jobject &quant_results,
		 size_t max_reads, bool plain_output) {

	// Initialize quantification class
	freq_orient_bias_quant_f fobquant (ref, alt);
	
	cerr << "ref: " << nuc_to_char(ref) << " alt: " << nuc_to_char(alt) << endl; 

	// Accumulate observed variables
	size_t n_reads = obquant_accumulate(fobquant, p, faidx, max_reads, plain_output);
	
	//calculate the estimators
	double theta_hat = fobquant.theta_hat();
	double phi_hat = fobquant();
	
	//output results to JSON
	using simpleson::json::jobject;
	jobject json_estimate;
	json_estimate["phi"] = phi_hat;

	jobject json_auxiliary;
	json_auxiliary["theta"] = theta_hat;

	jobject json_summary;
	json_summary["xc"] = fobquant.xc;
	json_summary["nc"] = fobquant.nc;
	json_summary["xi"] = fobquant.xi;
	json_summary["ni"] = fobquant.ni;

	quant_results["estimate"] = json_estimate;
	quant_results["auxiliary"] = json_auxiliary;
	quant_results["summary"] = json_summary;

	// output summary
	if (!plain_output) {
		frontend::global_log.v(2) << "xi: " << fobquant.xi << " xc: " << fobquant.xc << " ni: " << fobquant.ni <<
			" nc: " << fobquant.nc << '\n';
	}
	
	//output the estimator(s)
	if (!plain_output) {
		frontend::global_log.v(2) << "Theta estimator: ";
	}
	frontend::global_log.v(2) << theta_hat << '\n';
	if (!plain_output) {
		frontend::global_log.v(1) << "Phi estimator: ";
	} 
	frontend::global_log.v(1) << phi_hat << '\n';

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
bool orient_bias_quantify_bayes(nuc_t ref, nuc_t alt, piler &p, faidx_reader &faidx, simpleson::json::jobject &quant_results,
		 size_t max_reads, bool plain_output) {

	// Initialize objects
	bayes_orient_bias_quant_f bobquant (ref, alt);
	
	// accumulate observed variables
	size_t n_reads = obquant_accumulate(bobquant, p, faidx, max_reads, plain_output);
	
	// Get hyperparameter estimates (computationally expensive)
	hts::hparams result = bobquant.operator()<stograd::stepper::constant<double> > (BATCH_SIZE, N_EPOCHS, LRATE, CONVERGENCE_EPS);

	// output results to JSON
	using simpleson::json::jobject;
	jobject json_estimate;
	json_estimate["alpha_phi"] = result.alpha_phi;
	json_estimate["beta_phi"] = result.beta_phi;

	jobject json_auxiliary;
	json_auxiliary["alpha_theta"] = result.alpha_theta;
	json_auxiliary["beta_theta"] = result.beta_theta;

	jobject json_summary;
	json_summary["xc"] = accumulate(bobquant.m.xc_vec.begin(), bobquant.m.xc_vec.end(), 0);
	json_summary["nc"] = accumulate(bobquant.m.nc_vec.begin(), bobquant.m.nc_vec.end(), 0);
	json_summary["xi"] = accumulate(bobquant.m.xi_vec.begin(), bobquant.m.xi_vec.end(), 0);
	json_summary["ni"] = accumulate(bobquant.m.ni_vec.begin(), bobquant.m.ni_vec.end(), 0);

	quant_results += json_estimate;
	quant_results += json_auxiliary;
	quant_results += json_summary;

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

	// Output phi hyperparameters always, and E[phi] if verbosity >= 2
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