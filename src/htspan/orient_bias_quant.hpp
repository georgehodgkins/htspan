#ifndef _HTSPAN_ORIENT_BIAS_QUANT_HPP_
#define _HTSPAN_ORIENT_BIAS_QUANT_HPP_

#include <iostream>
#include <queue>
#include <vector>
#include <cassert>

#include <stdint.h>
#include <string.h>

#include <htslib/hts.h>

#include "bam.hpp"
#include "math.hpp"

namespace hts {

using namespace std;


struct orient_bias_quant_f {

	// oxoG-consistent alt count
	long int xc;
	
	// oxoG-inconsistent alt count
	long int xi;

	// oxoG-consistent total count
	long int nc;
	
	// oxoG-inconsistent total count
	long int ni;
	

	/// reference and alternative nucleotides to consider
	nuc_t r1_ref, r1_alt, r2_ref, r2_alt;

	/**
	 * Initialize class.
	 *
	 * Initialize class with expected ref and alt nucleotides for read 1
	 */
	orient_bias_quant_f(nuc_t _ref, nuc_t _alt)
	: r1_ref(_ref),
		r1_alt(_alt),
		r2_ref(nuc_complement(_ref)),
		r2_alt(nuc_complement(_alt))
	{
	}

	size_t size() const {
		return 0;
	}

	/**
	 * Push a read to accumulate statistics.
	 *
	 * @param b       BAM record of a read
	 * @param pos     reference position
	 */
	bool push(bam1_t* b, int32_t pos) {
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
			return false;
		}

		return true;
	}

	/**
	 * Push reads to accumulate statistics.
	 *
	 * @param pile bam pileup object
	 * @param n    number of reads in pileup object
	 * @param pos  target reference position
	 * @return number of successfully processed reads
	 */
	size_t push(const bam_pileup1_t* pile, size_t n, int32_t pos) {
		size_t success = 0;
		for (size_t i = 0; i < n; ++i) {
			if (push(pile[i].b, pos)) {
				++success;
			}
		}
		return success;
	}

	/**
	 * Return estimate of global DNA damage
	 *
	 * Reads should have been processed by calling `push`.
	 *
	 * @return estimate of phi
	 */
	double operator()() {
		double theta_hat = double(xi)/ni;
		return (double(xi)/nc - theta_hat) / (1 - theta_hat);
	}

};


}  // namespace hts

#endif  // _HTSPAN_ORIENT_BIAS_QUANT_HPP_
