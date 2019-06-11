#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <cstdlib>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "htspan/fetcher.hpp"
#include "htspan/orient_bias_filter.hpp"
#include "htspan/io/snv.hpp"

#include "htspan/io/simul_writer.hpp"
using namespace hts;

namespace hts {

/**
* This function provides an interface between the
* calling frontend and the backend classes for
* damage identification/filtering.
*
* @param ref Reference nucleotide of interest
* @param alt Alternate nucleotide of interest
* @param snv_fname Path to SNV file which locates SNVs to check for damage
* @param align_fname Path to BAM file containing sequence of interest
* @param phi Estimate of global damage rate
*/
bool orient_bias_filter(nuc_t ref, nuc_t alt, const char* snv_fname, const char* align_fname, double phi) {
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

			// run filter
			orient_bias_filter_f obfilter(ref, alt, f.pile.queries.size());
			obfilter.push(f.pile.queries, rec.pos, rec.nt_ref, rec.nt_alt);

			frontend::global_out << obfilter(phi) << '\t'
				<< obfilter.estimate_theta_given(phi, obfilter.estimate_initial_theta()) << '\n';
			f.clear();
		}
	}
	return true;
}

}// namespace hts