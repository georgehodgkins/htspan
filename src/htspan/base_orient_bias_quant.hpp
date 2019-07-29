#ifndef _HTSPAN_ORIENT_BIAS_QUANT_HPP_
#define _HTSPAN_ORIENT_BIAS_QUANT_HPP_

#include <vector>
#include <cassert>
#include <cstring>

#include <stdint.h>

#include <htslib/hts.h>

#include "bam.hpp"
#include "math.hpp"
#include "r_rand.hpp"
#include "nucleotide.hpp"

namespace hts {

using namespace std;

/**
* The base class for the frequentist and Bayesian quantification
* functors. Contains individual read push and simulation code.
*/
struct base_orient_bias_quant_f {

	/// reference and alternative nucleotides to consider
	nuc_t r1_ref, r1_alt, r2_ref, r2_alt;

	// oxoG-consistent alt count
	long int xc;
	
	// oxoG-inconsistent alt count
	long int xi;

	// oxoG-consistent total count
	long int nc;
	
	// oxoG-inconsistent total count
	long int ni;

	// count of total reads pushed
	int read_count;

	/**
	 * Initialize class.
	 *
	 * @param _ref Read 1 reference nucleotide for SNV of interest
	 * @param _alt Read 2 reference nucleotide for SNV of interest
	 * @param _p Initialized hts::piler containing sequence data of interest
	 * @param _f Initialized hts::faidx_reader containing reference sequence
	 */
	base_orient_bias_quant_f(nuc_t _ref, nuc_t _alt)
	: r1_ref(_ref),
		r1_alt(_alt),
		r2_ref(nuc_complement(_ref)),
		r2_alt(nuc_complement(_alt)),
		xc(0), xi(0), nc(0), ni(0), read_count(0)
	{
	}

	/*
	* The total number of reads pushed (not loci).
	* size() returns the number of loci instead.
	*/
	size_t n_reads() const {
		return read_count;
	}

	/**
	 * Process a single read and update observed variables accordingly.
	 *
	 * @param b       BAM record of a read
	 * @param pos     reference position
	 */
	bool push(const bam1_t* b, int32_t pos) {
		// only analyze nucleotides A, C, G, T (no indels)
		nuc_t qnuc = query_nucleotide(b, pos);
		if (!nuc_is_canonical(qnuc)) return false;
		// query aligning against the reverse strand is reverse-complemented,
		// so we need to reverse-complement again to get the original nucleotide
		nuc_t onuc;
		if (bam_is_rev(b)) {
			onuc = nuc_complement(qnuc);
		} else {
			onuc = qnuc;
		}

		// accmulate statistics based on read1 vs. read2 and original nucleotide
		if (bam_is_read1(b)) {
			if (onuc == r1_ref) {
				++nc;
			} else if (onuc == r1_alt) {
				++nc;
				++xc;
			} else if (onuc == r2_ref) {
				++ni;
			} else if (onuc == r2_alt) {
				++ni;
				++xi;
			}
			// other nucleotides are ignored
		// double-check that the read is second read, in case the flag is malformed
		} else if (bam_is_read2(b)) {
			if (onuc == r2_ref) {
				++nc;
			} else if (onuc == r2_alt) {
				++nc;
				++xc;
			} else if (onuc == r1_ref) {
				++ni;
			} else if (onuc == r1_alt) {
				++ni;
				++xi;
			}
			// other nucleotides are ignored
		} else {
			// read is neither read1 or read2
			// BAM file is likely not from paired end sequencing
			return false;
		}

		return true;
	}

	/*
	* Simulate a correctly distributed set of observed variables for a read, for testing.
	*
	* @param N Number of reads at the locus
	* @param theta Alternate allele probability at the locus
	* @param phi Damage probability at the locus
	* @param xc returned as damage-consistent alt count at the locus
	* @param xi returned as damage-inconsistent alt count at the locus
	* @param nc returned as correctly oriented alt count at the locus
	* @param ni returned as incorrectly oriented alt count at the locus
	*/
	static void simulate_orient_bias_read_counts (const size_t N, const double theta, const double phi,
			long int &xc, long int &xi, long int &nc, long int &ni, int seed = 0) {
		ni = N/2;
		nc = N - ni;
		xi = r_rand::rbinom(ni, theta, seed);
		long int xr = r_rand::rbinom(nc, theta, seed);
		long int xd = r_rand::rbinom(nc - xr, phi, seed);
		xc = xr + xd;
	}

	/**
	* Accumulate statistics for a set of reads at the given locus.
	*
	* @param pile Vector of pointers to the set of reads at the locus passing the quality filter
	* @param pos Reference position of the locus
	* @return Number of reads pushed.
	*/
	virtual size_t push (const vector<bam1_t*> &pile, int32_t pos) = 0;

	// access statistics for the last site pushed, for debug
	virtual long int xij () const = 0;
	virtual long int xcj () const = 0;
	virtual long int nij () const = 0;
	virtual long int ncj () const = 0;
};

}  // namespace hts

#endif  // _HTSPAN_ORIENT_BIAS_QUANT_HPP_
