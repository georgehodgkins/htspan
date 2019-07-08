#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <utility>
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
#include "htspan/io/simul_writer.hpp"
#include "htspan/io/streamer.hpp"

// TODO: update method documentation
namespace hts {

// Disambiguate non-fatal error codes from SNV reading (defined in snv.hpp)
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
			frontend::global_log.v(3) << "Info: fetched " << f.pile.queries.size() << " reads" << '\n';

			// read in data
			data.clear();
			data.reserve(f.pile.queries.size(), true);
			data.push(f.pile.queries, rec.pos, rec.nt_ref, rec.nt_alt);

		} else { // SNV is not damage-consistent
			frontend::global_log.v(3) << "Warning: Variant " << rec.to_string() <<
				" is not consistent with artifact type " << nuc_to_char(data.r1_ref) << '>' << nuc_to_char(data.r1_alt) <<
				'/' << nuc_to_char(data.r2_ref) << '>' << nuc_to_char(data.r2_alt) << ", ignoring.\n";
			snvr.err = -1;
		}
		return true;// not at EOF yet
	}
	return false;// no more records
}

/**
* This function provides an interface between the
* calling frontend and the backend classes for
* frequentist-model damage identification/filtering.
*
* @template SnvReader Derived class of snv::reader to use
* @param ref Reference nucleotide of interest
* @param alt Alternate nucleotide of interest
* @param eps Convergence tolerance for numerical minimization
* @param minz_bound Bounds for log-space numerical minimization
* @param snv_fname Path to SNV file which locates SNVs to check for damage
* @param align_fname Path to BAM file containing sequence of interest
* @param phi Estimate of global damage rate
* @return whether the operation succeeded.
*/
bool orient_bias_identify_freq(nuc_t ref, nuc_t alt, double eps, double minz_bound,
		 double phi, double sig_level, fetcher &f, snv::streamer &s, bool fixed_phi) {

	orient_bias_data data(ref, alt, 0);
	snv::record rec;
	snv::reader &snvr = *s.snvr_pt;
	snv::writer &snvw = *s.snvw_pt;

	// intialize filter
	// note that the filter is tied to the data object by reference, so it updates for new data
	freq_orient_bias_filter_f fobfilter(data, -minz_bound, minz_bound, eps);
	snvw.add_filter_tag(fobfilter);

	// table header
	frontend::global_log.v(1) << "snv\tpval\tfilter\n";
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
			}
			frontend::global_log.v(3) << " (skipped)\n";
			continue;
		}

		// run filter
		double pval = fobfilter(phi, fixed_phi);
		if (pval < sig_level) {
			++n_sig;
			frontend::global_log.v(2) << rec.to_string() << '\t' << pval << '\n';
			snvw.write(rec);
		} else {
			frontend::global_log.v(1) << rec.to_string() << '\t' << pval << "\t[fail]\n";
			snvw.write_filter_failed(rec, fobfilter);
		}
		++n;
	}
	frontend::global_log.v(1) << "\nTotal: " << n << " Passed: " << n_sig << " Failed: " << n - n_sig << '\n';
	return true;
}

/**
* This function provides an interface between the
* calling frontend and the backend classes for
* Bayesian-model damage identification/filtering.
*
* @template SnvReader Derived class of snv::reader to use
* @param ref Reference nucleotide of interest
* @param alt Alternate nucleotide of interest
* @param eps Convergence tolerance for numerical minimization
* @param minz_bound Bounds for log-space numerical minimization
* @param snv_fname Path to SNV file which locates SNVs to check for damage
* @param align_fname Path to BAM file containing sequence of interest
* @param alpha Alpha hyperparameter for the beta distribution
* @param beta Beta hyperparameter for the beta distribution
* @param prior_alt Prior probability of the alternate model (theta != 0)
* @return whether the operation succeeded.
*/
bool orient_bias_identify_bayes(nuc_t ref, nuc_t alt, double eps, double minz_bound, 
		double alpha, double beta, double prior_alt, double posterior_threshold, fetcher &f, snv::streamer &s) {
	orient_bias_data data (ref, alt, 0);
	snv::record rec;
	snv::reader &snvr = *s.snvr_pt;
	snv::writer &snvw = *s.snvw_pt;

	// intialize filter
	// note that the filter is tied to the data object by reference, so it updates for new data
	bayes_orient_bias_filter_f bobfilter(data, -minz_bound, minz_bound, eps);
	snvw.add_filter_tag(bobfilter);

	// table header
	frontend::global_log.v(1) << "snv\tprob\tfilter" << '\n';
	size_t n_sig = 0;
	size_t n = 0;

	// note that fetch_next_snv modifies all of the objects passed to it
	// In particular, data is populated with a new set of data for the read SNV
	while (fetch_next_snv(snvr, f, data, rec)) {
			// check for non-fatal errors in SNV reading/piling
			if (snvr.error()) {
				if (snvr.error() == -1) {// read suceeded, inconsistent variants
					snvw.write(rec);
				}
				continue;
			}

			// run filter
			double lposterior = bobfilter(prior_alt, alpha, beta);
			double posterior = exp(lposterior);
			if (posterior > posterior_threshold) {
				++n_sig;
				frontend::global_log.v(2) << rec.to_string() << '\t' << posterior << '\n';
				snvw.write(rec);
			} else {
				frontend::global_log.v(1) << rec.to_string() << '\t' << posterior << "\t[fail]\n";
				snvw.write_filter_failed(rec, bobfilter);
			}
			++n;
	}

	frontend::global_log.v(1) << "\nSummary: Total: " << n << " Passed: " << n_sig << " Failed: " << n - n_sig << '\n';

	return true;
}

}// namespace hts