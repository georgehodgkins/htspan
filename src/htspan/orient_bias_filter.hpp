#ifndef _HTSPAN_ORIENT_BIAS_HPP_
#define _HTSPAN_ORIENT_BIAS_HPP_

#include <iostream>
#include <queue>
#include <vector>
#include <cassert>

#include <stdint.h>
#include <string.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "fetcher.hpp"
#include "math.hpp"

namespace hts {

using namespace std;

struct orient_bias_filter_f {

	/// reference base (0), alternative base (1), or other base (2)
	vector<int8_t> bases;
	/// base error probability
	vector<double> errors;
	/// oxidation probability: phi if orientation is consistent, 0 otherwise
	vector<double> phis;

	/// global oxidation error probablity
	double phi;

	/// reference and alternative nucleotides to consider
	nucleotide f_ref, f_alt, r_ref, r_alt;

	/// observed nucleotide of the read
	vector<char> cnucs;
	/// base quality score at the query position
	vector<int> quals;
	/// forward or reverse reference strand to which the read aligned
	vector<char> strands;
	/// first read or second read of a pair
	vector<int> members;

	orient_bias_filter_f(nucleotide _ref, nucleotide _alt, double _phi, size_t n)
	: f_ref(_ref),
		f_alt(_alt),
		r_ref((nucleotide)nuc_complement(_ref)),
		r_alt((nucleotide)nuc_complement(_alt)),
		phi(_phi)
	{
		bases.reserve(n);
		errors.reserve(n);
		phis.reserve(n);
	}

	size_t size() const {
		return bases.size();
	}

	/**
	 * Push a read to accumulate statistics.
	 */
	bool push(bam1_t* b, int32_t pos) {
		nucleotide qnuc = (nucleotide) query_nucleotide(b, pos);
		// only analyze nucleotides A, C, G, T (no indels)
		if (nuc_is_canonical(qnuc)) {
			// get the original nucleotide of the read
			if (bam_is_rev(b)) {
				qnuc = (nucleotide) nuc_complement(qnuc);
			}
			if (bam_is_read1(b)) {
				if (qnuc == f_ref) {
					// e.g. G>G on read 1
					bases.push_back(0);
					phis.push_back(0);
				} else if (qnuc == f_alt) {
					// e.g. G>T on read 1
					bases.push_back(1);
					phis.push_back(phi);
				} else {
					// e.g. G>A on read 1
					// query nucleotide is a non-ref, non-alt nucleotide
					bases.push_back(2);
					phis.push_back(0);
				}
			// double-check that the read is second read, in case the flag is malformed
			} else if (bam_is_read2(b)) {
				if (qnuc == r_ref) {
					// e.g. C>C on read 2
					bases.push_back(0);
					phis.push_back(0);
				} else if (qnuc == r_alt) {
					// e.g. C>A on read 2
					bases.push_back(1);
					phis.push_back(phi);
				} else {
					// e.g. C>G or C>T on read 2
					bases.push_back(2);
					phis.push_back(0);
				}
			}  // if (bam_is_read1(b))

			errors.push_back(anti_phred((double)query_quality(b, pos)));

			// populate informational fields
			
			quals.push_back(query_quality(b, pos));

			if (bam_is_rev(b)) {
				strands.push_back('-');
			} else {
				strands.push_back('+');
			}

			if (bam_is_read1(b)) {
				members.push_back(1);
			} else if (bam_is_read2(b)) {
				members.push_back(2);
			} else {
				members.push_back(0);
			}

			switch (qnuc) {
			case nuc_A:
				cnucs.push_back('A');
				break;
			case nuc_C:
				cnucs.push_back('C');
				break;
			case nuc_G:
				cnucs.push_back('G');
				break;
			case nuc_T:
				cnucs.push_back('T');
				break;
			}

		}  // if (nuc_is_canonical(qnuc))
	}

};

}  // namespace hts

#endif  // _HTSPAN_ORIENT_BIAS_HPP_
