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
#include "htspan/io/snv.hpp"

#include "htspan/io/simul_writer.hpp"
using namespace hts;

namespace hts {

/**
* This function provides an interface between the
* calling frontend and the backend classes for
* frequentist-model damage identification/filtering.
*
* @param ref Reference nucleotide of interest
* @param alt Alternate nucleotide of interest
* @param snv_fname Path to SNV file which locates SNVs to check for damage
* @param align_fname Path to BAM file containing sequence of interest
* @param phi Estimate of global damage rate
* @return whether the operation succeeded.
*/
bool orient_bias_identify_freq(nuc_t ref, nuc_t alt, const char* snv_fname, const char* align_fname, double phi) {
	fetcher f;
	snv::reader snvr (snv_fname); // throws an exception if opening fails
	snv::record rec;
	// Open BAM file
	if (!f.open(align_fname)) {
		std::cerr <<
			"Error: Could not open BAM file \'" << align_fname << "\'.";
		return false;
	}
	// table header
	frontend::global_out << "p\ttheta_hat" << '\n';
	while (snvr.next(rec)) {
		//check that the SNV is consistent with damage
		if ((rec.nt_ref == ref && rec.nt_alt == alt) ||
				(rec.nt_ref == nuc_complement(ref) && rec.nt_alt == nuc_complement(alt))) {
			// fetch reads at target position
			int32_t rid = bam_name_to_id(f.hdr, rec.chrom);
			if (rid == -1) {
				frontend::global_log.v(1) << "Warning: could not find " << rec.chrom << '\n';
				continue;
			}

			if (!f.fetch(rid, rec.pos)) {
				frontend::global_log.v(1) << "Warning: could not fetch reads for: " << rec.chrom << ':' << rec.pos << '\n';
				continue;
			}

			frontend::global_log.v(2) << "Info: fetched " << f.pile.queries.size() << " reads" << '\n';

			// read in data
			orient_bias_data data(ref, alt, f.pile.queries.size());
			data.push(f.pile.queries, rec.pos, rec.nt_ref, rec.nt_alt);

			// run filter
			freq_orient_bias_filter_f fobfilter(data);
			frontend::global_out << fobfilter(phi) << '\t'
				<< fobfilter.estimate_theta_given(phi, fobfilter.estimate_initial_theta()) << '\n';
			f.clear();
		}
	}
	return true;
}

/**
* This function provides an interface between the
* calling frontend and the backend classes for
* Bayesian-model damage identification/filtering.
*
* @template Integration method to use
* @param ref Reference nucleotide of interest
* @param alt Alternate nucleotide of interest
* @param snv_fname Path to SNV file which locates SNVs to check for damage
* @param align_fname Path to BAM file containing sequence of interest
* @param phi Estimate of global damage rate
* @return whether the operation succeeded.
*/
template <typename Integrator>
bool orient_bias_identify_bayes(nuc_t ref, nuc_t alt, const char* snv_fname, const char* align_fname,
		double alpha, double beta, double prior_alt) {
	fetcher f;
	snv::reader snvr (snv_fname); // throws an exception if opening fails
	snv::record rec;
	// Open BAM file
	if (!f.open(align_fname)) {
		std::cerr <<
			"Error: Could not open BAM file \'" << align_fname << "\'.";
		return false;
	}
	// table header
	frontend::global_out << "log posterior prob. of non-zero theta" << '\n';
	while (snvr.next(rec)) {
		//check that the SNV is consistent with damage
		if ((rec.nt_ref == ref && rec.nt_alt == alt) ||
				(rec.nt_ref == nuc_complement(ref) && rec.nt_alt == nuc_complement(alt))) {
			// fetch reads at target position
			int32_t rid = bam_name_to_id(f.hdr, rec.chrom);
			if (rid == -1) {
				frontend::global_log.v(1) << "Warning: could not find " << rec.chrom << '\n';
				continue;
			}

			if (!f.fetch(rid, rec.pos)) {
				frontend::global_log.v(1) << "Warning: could not fetch reads for: " << rec.chrom << ':' << rec.pos << '\n';
				continue;
			}

			frontend::global_log.v(2) << "Info: fetched " << f.pile.queries.size() << " reads" << '\n';

			// run filter
			orient_bias_data data(ref, alt, f.pile.queries.size());
			data.push(f.pile.queries, rec.pos, rec.nt_ref, rec.nt_alt);

			// run filter
			bayes_orient_bias_filter_f bobfilter(data);
			frontend::global_out << bobfilter.operator()<Integrator>(prior_alt, alpha, beta);
			f.clear();
		}
	}
	return true;
}

}// namespace hts