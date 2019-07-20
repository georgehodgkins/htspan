#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <limits>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "htspan/fetcher.hpp"
#include "htspan/freq_orient_bias_filter.hpp"
#include "htspan/bayes_orient_bias_filter.hpp"
#include "htspan/nucleotide.hpp"

#include "htspan/io/snv_reader.hpp"
#include "htspan/io/snv_writer.hpp"
#include "htspan/io/snv.hpp"
#include "htspan/io/streamer.hpp"

#include "simul_writer.hpp"

/**
* This file contains driver functions that implement the identification process.
*
* Alignment data from a hts::fetcher and SNV data from a hts::snv::snvr are collated and
* stored in an orient_bias_data object, which is then (for damage-consistent variants) linked to
* a filter class (Bayesian or frequentist), which returns a statistic indicating how likely it is
* that each variant is the result of damage. These SNVs are written back to an SNV file, with a field
* recording the statistic and optionally, for VCF files, an entry in the filter column if they failed.
*/ 

namespace hts {

// Disambiguate non-fatal error codes from SNV reading (defined in snv.hpp)
// TODO move to snv.hpp
inline void print_snvr_err(snv::reader &snvr) {
	switch(snvr.error()) {
	case -2:
		frontend::global_log.v(1) << "Warning: Reader has not been initialized.\n";
		return;
	case 1:
		frontend::global_log.v(1) << "Warning: ignoring zero-nucleotide/invalid variant read.\n";
		return;
	case 2:
		frontend::global_log.v(1) << "Warning: could not unpack VCF/BCF record.\n";
		return;
	case 3:
		frontend::global_log.v(1) << "Warning: ignoring multiple-nucleotide variant read.\n";
		return;
	case 0:
	default:
		return;
	}
}

/**
* This function reads the next SNV from the given SNV reader,
* fetches corresponding reads from the given fetcher, and piles
* the fetched reads in the given data object.
*
* @param snvr Any derived class of snv::reader, already opened
* @param f A hts::fetcher, already opened
* @param data orient_bias_data object initialized with the nucleotides of interest
* @param rec snv::record to read data into
* @return Whether the snvr has reached EOF
*	All passed objects are also modified:
*		snvr & f internal iterators are advanced
*		data is cleared and populated with reads corresponding to the read SNV
*		rec is populated with the read SNV, including error codes if reading or fetching fails
*/
bool fetch_next_snv (snv::reader &snvr, fetcher &f, orient_bias_data &data, snv::record &rec) {
	if (snvr.next(rec)) {
		
		// check for non-fatal errors
		if (snvr.error() != 0) {
			print_snvr_err(snvr);
			return true;
		}
		
		// check that the SNV is consistent with damage
		if ((rec.nt_ref == data.r1_ref && rec.nt_alt == data.r1_alt) ||
				(rec.nt_ref == data.r2_ref && rec.nt_alt == data.r2_alt)) {
					
			// get ref ID from chromosome string
			int rid = bam_name_to_id(f.hdr, rec.chrom);
			if (rid == -1) {
				frontend::global_log.v(1) << "Warning: could not find " << rec.chrom << '\n';
				return true;
			}
			
			// fetch reads at target position
			f.clear();
			if (!f.fetch(rid, rec.pos)) {
				frontend::global_log.v(1) << "Warning: could not fetch reads for " << rec.to_string() << '\n';
				snvr.err = -1;
				return true;
			}

			// read in data
			data.clear();
			data.reserve(f.pile.queries.size(), true);
			data.push(f.pile.queries, rec.pos, rec.nt_ref, rec.nt_alt);

		} else { // SNV is not damage-consistent
			snvr.err = -1;
		}
		
		return true; // not at EOF yet
	}
	
	return false; // no more records
}

/**
* This function provides a driver for the process of frequentist identification.
*
* Verbosity output levels:
* 	1: warnings and failing (non-significant) SNVs and p-vals,
*		or all damage-consistent SNVS if significance is not considered
* 	2: number of reads fetched per SNV, significant (passing) SNVs if significance is being considered
* 	3: Damage-inconsistent SNVs
*
*
* @param ref Reference nucleotide of interest
* @param alt Alternate nucleotide of interest
* @param eps Convergence tolerance for numerical minimization
* @param minz_bound Bounds for log-space numerical minimization
* @param phi Initial estimate of global damage rate
* @param sig_level Maximum p-val considered significant (i.e. variant is genuine);
*	0 indicates significance not considered
* @param f Initialized fetcher containing the alignment of interest
* @param s Initialized streamer containing input and output SNV files
* @param fixed_phi Whether to try to improve our estimate of phi during identification
* @param plain_output Whether to output non machine-readable portions of output
* @return whether the operation succeeded.
*/

bool orient_bias_identify_freq(nuc_t ref, nuc_t alt, double eps, double minz_bound,
		 double phi, double sig_level, fetcher &f, snv::streamer &s, bool fixed_phi, bool plain_output) {

	orient_bias_data data(ref, alt, 0);
	snv::record rec;
	snv::reader &snvr = *s.snvr_pt;
	snv::writer &snvw = *s.snvw_pt;

	// add a data field for p-values (column in TSV, INFO field in VCF)
	snvw.add_numeric_info("FOBP", "P-value for genuine variant from hts-orient-bias filter");

	// intialize filter
	// note that the filter is tied to the data object by reference, so it updates for new data
	freq_orient_bias_filter_f fobfilter(data, -minz_bound, minz_bound, eps);
	if (sig_level > 0) {
		snvw.add_filter_tag(fobfilter);
	}

	// table header
	frontend::global_log.v(1) << "snv\tpval";
	if (sig_level > 0) {
		frontend::global_log.v(1) << "\tfilter";
	}
	frontend::global_log.v(2) << "\t#reads";
	frontend::global_log.v(1) << "\n";
	size_t n_sig = 0;
	size_t n = 0;

	// note that fetch_next_snv modifies all of the objects passed to it
	// In particular, data is populated with a new set of data for the read SNV
	while (fetch_next_snv(snvr, f, data, rec)) {
		// check for non-fatal errors in SNV reading/piling and
		// skip this record if an error was encountered
		// (fetch_next_snv should print the appropriate warning)
		if (snvr.error()) {
			if (snvr.error() == -1) { // successful read, inconsistent variant
				snvw.write(rec);
				frontend::global_log.v(3) << rec.to_string() << "\tN/A\t[nc]\tN/A\n";
			}
			continue;
		}

		// run filter
		double pval = fobfilter(phi, fixed_phi);
		
		if (sig_level == 0) { // significance is not considered
			frontend::global_log.v(1) << rec.to_string() << '\t' << pval;
			frontend::global_log.v(2) << '\t' <<  f.pile.queries.size();
			frontend::global_log.v(1) << '\n';
			snvw.write(rec, "FOBP", pval);
			
		} else if (pval < sig_level) { // variant is significant (probably genuine)
			++n_sig;
			frontend::global_log.v(2) << rec.to_string() << '\t' << pval << "\t[pass]\t" << f.pile.queries.size() << '\n';
			snvw.write(rec, "FOBP", pval);
			
		} else { // variant is not significant (probably damage)
			frontend::global_log.v(1) << rec.to_string() << '\t' << pval << "\t[fail]";
			frontend::global_log.v(2) << '\t' <<  f.pile.queries.size();
			frontend::global_log.v(1) << '\n';
			snvw.write_filter_failed(rec, fobfilter, "FOBP", pval);
		}
		
		++n;
	}
	
	// Extra caveats and summary data that are not important to a computer
	if (!plain_output) {
		frontend::global_log.v(1) << "NOTE 1: P-values are not adjusted for false discovery.\n";
		frontend::global_log.v(1) << "NOTE 2: Statistic values are only accurate to within " << eps << ",\n\t" <<
			"so p-vals of 0 or 1 are not true 0 or 1";
			
		frontend::global_log.v(1) << "\nTotal: " << n;
		
		if (sig_level > 0) {
			frontend::global_log.v(1) << " Passed: " << n_sig << " Failed: " << n - n_sig;
		}
		
		frontend::global_log.v(1) << '\n';
	}
	return true;
}

/**
* This function provides a driver for the Bayesian identification process
*
* Verbosity output levels:
* 	1: warnings and failing (non-significant) SNVs and posterior probs,
*		or all damage-consistent SNVS if significance is not considered
* 	2: number of reads fetched per SNV, significant (passing) SNVs if significance is being considered
* 	3: Damage-inconsistent SNVs
*
* @template SnvReader Derived class of snv::reader to use
* @param ref Reference nucleotide of interest
* @param alt Alternate nucleotide of interest
* @param eps Convergence tolerance for numerical minimization
* @param minz_bound Bounds for log-space numerical minimization
* @param alpha Alpha hyperparameter for the beta distribution
* @param beta Beta hyperparameter for the beta distribution
* @param prior_alt Prior probability of the alternate model (theta != 0)
* @param posterior_threshold Posterior probability of alt model must be above this threshold
*	for a variant to be considered significant; 0 indicates significance is not considered
* @param f Already initialized fetcher containing the alignment of interest
* @param s Already initialized streamer containing input and output SNV files
* @param plain_output Whether to print non machine-readable portions of output
* @return whether the operation succeeded.
*/
bool orient_bias_identify_bayes(nuc_t ref, nuc_t alt, double eps, double minz_bound, 
		double alpha, double beta, double prior_alt, double posterior_threshold, fetcher &f, snv::streamer &s, bool plain_output) {
	orient_bias_data data (ref, alt, 0);
	snv::record rec;
	snv::reader &snvr = *s.snvr_pt;
	snv::writer &snvw = *s.snvw_pt;

	// intialize filter
	// note that the filter is tied to the data object by reference, so it updates for new data
	bayes_orient_bias_filter_f bobfilter(data, -minz_bound, minz_bound, eps);
	if (posterior_threshold > 0) {
		snvw.add_filter_tag(bobfilter);
	}

	// add info field for posterior probabilities
	snvw.add_numeric_info("BOBP", "Log posterior probability of a genuine variant from hts-orient-bias filter");

	// table header
	frontend::global_log.v(1) << "snv\tprob";
	if (posterior_threshold > 0) {
		frontend::global_log.v(1) << "\tfilter";
	}
	frontend::global_log.v(2) << "\t#reads";
	frontend::global_log.v(1) << '\n';
	
	size_t n_sig = 0; // number of significant variants
	size_t n = 0; // number of analyzed variants

	// note that fetch_next_snv modifies all of the objects passed to it
	// In particular, data is populated with a new set of data for the read SNV
	while (fetch_next_snv(snvr, f, data, rec)) {
			// check for non-fatal errors in SNV reading/piling
			if (snvr.error()) {
				if (snvr.error() == -1) {// read suceeded, inconsistent variants
					frontend::global_log.v(3) << rec.to_string() << "\tN/A\t[nc]\tN/A\n";
					snvw.write(rec);
				}
				continue;
			}
			
			// run filter (returns log of posterior prob of alt model)
			double lposterior = bobfilter(prior_alt, alpha, beta);
			double posterior = exp(lposterior);

			if (posterior_threshold == 0) { // significance is not considered
				frontend::global_log.v(1) << rec.to_string() << '\t' << posterior;
				frontend::global_log.v(2) << '\t' <<  f.pile.queries.size();
				frontend::global_log.v(1) << '\n';
				snvw.write(rec, "BOBP", lposterior);
				
			} else if (posterior > posterior_threshold) { // variant is significant (likely genuine)
				++n_sig;
				frontend::global_log.v(2) << rec.to_string() << '\t' << posterior << "\t[pass]\t" << f.pile.queries.size() << '\n';
				snvw.write(rec, "BOBP", lposterior);
				
			} else { // variant is not significant (likely damage)
				frontend::global_log.v(1) << rec.to_string() << '\t' << posterior << "\t[fail]";
				frontend::global_log.v(2) << '\t' <<  f.pile.queries.size();
				frontend::global_log.v(1) << '\n';
				snvw.write_filter_failed(rec, bobfilter, "BOBP", lposterior);
			}
			++n; // number of analyzed variants
	}
	
	// Additional caveats and other non machine-readable output
	if (!plain_output) {
		
		frontend::global_log.v(1) << "NOTE 1: posterior probabilities are not adjusted for false discovery.\n";
		frontend::global_log.v(1) << "NOTE 2: Values are only accurate to within " << eps << ",\n\t" <<
			"so probabilities of 0 or 1 are not true 0 or 1.";

		frontend::global_log.v(1) << "\nSummary: Total: " << n;
		if (posterior_threshold > 0) {
			frontend::global_log.v(1) << " Passed: " << n_sig << " Failed: " << n - n_sig;
		}
		frontend::global_log.v(1) << '\n';
	}

	return true;
}

}// namespace hts