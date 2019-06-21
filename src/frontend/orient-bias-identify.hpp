#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <cstdlib>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "htspan/fetcher.hpp"
#include "htspan/freq_orient_bias_filter.hpp"
#include "htspan/bayes_orient_bias_filter.hpp"
#include "htspan/nucleotide.hpp"

#include "htspan/io/snv_reader.hpp"
#include "htspan/io/simul_writer.hpp"


namespace hts {

// Disambiguate non-fatal error codes from SNV reading (defined in snv.hpp)
inline void print_snvr_err(snv::record &rec) {
	switch(rec.err) {
	case 1:
		frontend::global_log.v(1) << "Warning: could not find " << rec.chrom;
		return;
	case 2:
		frontend::global_log.v(1) << "Warning: could not unpack VCF/BCF record.";
		return;
	case 3:
		frontend::global_log.v(1) << "Warning: ignoring multiple-nucleotide variant read from VCF/BCF.";
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
		if (rec.err != 0) {
			print_snvr_err(rec);
			return true;
		}
		// check that the SNV is consistent with damage
		if ((rec.nt_ref == data.r1_ref && rec.nt_alt == data.r1_alt) ||
				(rec.nt_ref == data.r2_ref && rec.nt_alt == data.r2_alt)) {
			// get ref ID from chromosome string
			int rid = bam_name_to_id(f.hdr, rec.chrom);
			if (rid == -1) {
				frontend::global_log.v(1) << "Warning: could not find " << rec.chrom;
				return true;
			}
			// fetch reads at target position
			f.clear();
			if (!f.fetch(rid, rec.pos)) {
				frontend::global_log.v(1) << "Warning: could not fetch reads for: " << rec.chrom << ':' << rec.pos << '\n';
				rec.err = -1;
				return true;
			}
			frontend::global_log.v(2) << "Info: fetched " << f.pile.queries.size() << " reads" << '\n';

			// read in data
			data.clear();
			data.reserve(f.pile.queries.size(), true);
			data.push(f.pile.queries, rec.pos, rec.nt_ref, rec.nt_alt);

		} else { // SNV is not damage-consistent
			frontend::global_log.v(1) << "Warning: Variant " << nuc_to_char(rec.nt_ref) << '>' <<
				nuc_to_char(rec.nt_alt) << "is not consistent with selected damage type, ignoring.";
			rec.err = -1;
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
template <typename SnvReader>
bool orient_bias_identify_freq(nuc_t ref, nuc_t alt, double eps, double minz_bound,
		 const char* snv_fname, const char* align_fname, double phi) {
	fetcher f;
	orient_bias_data data(ref, alt, 0);
	// Open BAM file
	if (!f.open(align_fname)) {
		std::cerr <<
			"Error: Could not open BAM file \'" << align_fname << "\'.";
		return false;
	}
	// open SNV file (TSV or VCF)
	SnvReader snvr_obj(snv_fname); // throws an exception if opening fails
	snv::reader &snvr = snvr_obj; // assign a base class reference to the derived class object

	snv::record rec;
	// table header
	frontend::global_out << "p\ttheta_hat" << '\n';
	// note that pile_next_snv modifies all of the objects passed to it
	// In particular, data is populated with a new set of data for the read SNV
	while (fetch_next_snv(snvr, f, data, rec)) {
			// check for non-fatal errors in SNV reading/piling
			if (rec.err) {
				continue;
			}
			// run filter
			freq_orient_bias_filter_f fobfilter(data, -minz_bound, minz_bound, eps);
			frontend::global_out << fobfilter(phi) << '\t'
				<< fobfilter.estimate_theta_given(phi) << '\n';
		}
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
template <typename SnvReader>
bool orient_bias_identify_bayes(nuc_t ref, nuc_t alt, double eps, double minz_bound, const char* snv_fname,
		const char* align_fname, double alpha, double beta, double prior_alt) {
	fetcher f;
	orient_bias_data data (ref, alt, 0);
	// Open BAM file
	if (!f.open(align_fname)) {
		std::cerr <<
			"Error: Could not open BAM file \'" << align_fname << "\'.";
		return false;
	}
	// open SNV file (TSV or VCF)
	SnvReader snvr_obj(snv_fname);// throws exception on failure to open
	snv::reader &snvr = snvr_obj;
	snv::record rec;
	// table header
	frontend::global_out << "log posterior prob. of non-zero theta" << '\n';
	// note that pile_next_snv modifies all of the objects passed to it
	// In particular, data is populated with a new set of data for the read SNV
	while (fetch_next_snv(snvr, f, data, rec)) {
			// check for non-fatal errors in SNV reading/piling
			if (rec.err) {
				continue;
			}
			// run filter
			bayes_orient_bias_filter_f bobfilter(data, -minz_bound, minz_bound, eps);
			frontend::global_out << bobfilter(prior_alt, alpha, beta);
	}
	return true;
}

}// namespace hts