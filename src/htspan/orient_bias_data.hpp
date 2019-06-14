#ifndef _HTSPAN_ORIENT_BIAS_DATA_HPP_
#define _HTSPAN_ORIENT_BIAS_DATA_HPP_

#include <vector>
#include <stdint.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <gsl/gsl_randist.h>// for simulation

#include "io/orient_bias_stats.hpp"
#include "fetcher.hpp"
#include "nucleotide.hpp"

namespace hts {

struct orient_bias_data {

	/// reference base (0), alternative base (1), or other base (2)
	vector<int8_t> bases;
	/// base error probability
	vector<double> errors;
	/// whether variant has an orientation consistent with artifact damage
	/// (false if base is not variant)
	vector<bool> orients;

	/// observed nucleotide of the read
	vector<char> cnucs;
	/// base quality score at the query position
	vector<int> quals;
	/// forward or reverse reference strand to which the read aligned
	vector<char> strands;
	/// first read or second read of a pair
	vector<int> members;

	/// reference and alternative nucleotides to consider
	nuc_t r1_ref, r1_alt, r2_ref, r2_alt;


	orient_bias_data (nuc_t _ref, nuc_t _alt, size_t n)
		: r1_ref(_ref),
			r1_alt(_alt),
			r2_ref(nuc_complement(r1_ref)),
			r2_alt(nuc_complement(r2_alt))
	{
		reserve(n);
	}

	/**
	 * Push a read to accumulate statistics.
	 *
	 * @param b       BAM record of a read
	 * @param pos     reference position
	 * @param nt_ref  reference nucleotide
	 * @param nt_alt  alternative nucleotide
	 */
	bool push(bam1_t* b, int32_t pos, nuc_t nt_ref, nuc_t nt_alt) {
		// only analyze nucleotides A, C, G, T (no indels)
		nuc_t qnuc = query_nucleotide(b, pos);
		if (!nuc_is_canonical(qnuc)) return false;

		bool success = true;

		// determine if query nucleotide matches reference or alternative
		// NB  query aligning against the reverse strand is reverse-complemented
		//     s.t. its sequence is in the same direction as the reference;
		//     the same is true for the reported nt_ref and nt_alt
		if (qnuc == nt_ref) {
			bases.push_back(0);
		} else if (qnuc == nt_alt) {
			bases.push_back(1);
		} else {
			bases.push_back(2);
		}

		// call orientation (damage-consistency) based on the original nucleotide of the read
		nuc_t onuc;
		if (bam_is_rev(b)) {
			onuc = nuc_complement(qnuc);
		} else {
			onuc = qnuc;
		}
		if (bam_is_read1(b)) {
			if (onuc == r1_ref || onuc == r1_alt) {
				// e.g. G or T on read 1
				orients.push_back(true);
			} else {
				orients.push_back(false);
			}
		// double-check that the read is second read, in case the flag is malformed
		} else if (bam_is_read2(b)) {
			members.push_back(2);
			if (onuc == r2_ref || onuc == r2_alt) {
				// e.g. C or A on read 2
				orients.push_back(true);
			} else {
				orients.push_back(false);
			}
		} else {
			// flag is malformed
			members.push_back(0);
			orients.push_back(false);
			success = false;
		}

		errors.push_back(anti_phred((double)query_quality(b, pos)));

		// populate informational fields
		
		quals.push_back(query_quality(b, pos));

		if (bam_is_rev(b)) {
			strands.push_back('-');
		} else {
			strands.push_back('+');
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
		default:
			// this should never happen
			cnucs.push_back('N');
			break;
		}

		if (!success) {
			// rollback everything
			bases.pop_back();
			errors.pop_back();
			orients.pop_back();
			cnucs.pop_back();
			quals.pop_back();
			strands.pop_back();
			members.pop_back();
		}

		return success;
	}

	/**
	 * Push reads to accumulate statistics.
	 *
	 * @param bs   pile of reads
	 * @param pos  target reference position
	 * @param nt_ref  reference nucleotide
	 * @param nt_alt  alternative nucleotide
	 * @return number of successfully processed reads
	 */
	size_t push(vector<bam1_t*> bs, int32_t pos, nuc_t nt_ref, nuc_t nt_alt) {
		size_t success = 0;
		for (size_t i = 0; i < bs.size(); ++i) {
			if (push(bs[i], pos, nt_ref, nt_alt)) {
				++success;
			}
		}
		return success;
	}

	/**
	* C++ version of the read simulation code in R,
	* populating class members without real BAM data.
	* Used to make unit tests consistent with R model.
	*
	* Note that this will clear any reads previously pushed.
	*
	* @param theta Target theta value
	* @param phi Target phi value
	* @param err_mean Mean value for normal distribution of error values
	* @param err_sd Std dev for normal distribution of error values
	* @return none (populates class members with simulated values)
	*/
	void simulate(double theta, double phi, double err_mean, double err_sd) {
		size_t count = bases.capacity();//set in constructor
		clear();
		reserve(count);//in case clearing resets allocation, which is implementation dependent
		//NB: R docs are unclear about what RNG is used for random sampling, so this may be a point of divergence
		gsl_rng *rng = gsl_rng_alloc (gsl_rng_mt19937);
		//generate normally distributed error values with the input parameters
		for (size_t n = 0; n < count; ++n) {
			double z = gsl_ran_gaussian(rng, err_sd) + err_mean;
			if (z < 0) z = 0;
			errors[n] = anti_phred(z);
		}
		//generate bases where prob of alt is theta and prob of ref is 1-theta
		for (size_t n = 0; n < count; ++n) {
			double z = gsl_ran_flat(rng, 0, 1);
			bases[n] = (z < theta) ? 0 : 1;//0=ref, 1=alt
			//generate orients with equal probability
			z = gsl_ran_flat(rng, 0, 1);
			orients[n] = (z < 0.5);
			//for oriented reference bases, prob of damage (ref>alt) is phi
			if (bases[n] == 0 && orients[n]) {
				z = gsl_ran_flat(rng, 0, 1);
				if (z < phi) {
					bases[n] = 1;
				}
			}
		}
	}

	/**
	* Reserve space in vector members for a given number of reads.
	*
	* @param n Number of fields to reserve in each member
	* @param info_fields Whether to reserve space in informational vectors as well [false]
	*/
	void reserve (size_t n, bool info_fields = false) {
		orients.reserve(n);
		bases.reserve(n);
		errors.reserve(n);
		if (info_fields) {
			cnucs.reserve(n);
			quals.reserve(n);
			strands.reserve(n);
			members.reserve(n);
		}
	}

	/**
	* Clear member vectors.
	*/
	void clear () {
		orients.clear();
		bases.clear();
		errors.clear();
		cnucs.clear();
		quals.clear();
		strands.clear();
		members.clear();
	}

	/**
	* Reads previously generated simulation data from a .tsv file,
	* and populates class members with that data.
	*
	* Note that this will clear any reads previously pushed.
	*
	* @param path to the .tsv file
	* @return none (populates orients, bases, and errors)
	*/
	void read (const char* path) {
		clear();
		orient_bias_stats::reader stat_reader (path);
		orient_bias_stats::record rec;
		while (stat_reader.next(rec)) {
			bases.push_back(rec.base);
			errors.push_back(rec.error);
			orients.push_back(rec.orient);
		}
		stat_reader.close();
	}
};

}// namespace hts

#endif // _HTSPAN_ORIENT_BIAS_DATA_HPP_ 