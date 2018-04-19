#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <cstdlib>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "htspan/fetcher.hpp"
#include "htspan/math.hpp"
using namespace hts;

/**
 * Get whether query is read 1.
 * 
 * @param b  pointer to a BAM record
 * @return true if query is read 1
 */
#define bam_is_read1(b) (((b)->core.flag & BAM_FREAD1) != 0)

/**
 * Get whether query is read 2.
 * 
 * @param b  pointer to a BAM record
 * @return true if query is read 1
 */
#define bam_is_read2(b) (((b)->core.flag & BAM_FREAD2) != 0)

int main(int argc, char** argv) {
	--argc;

	if (argc != 7) {
		cerr << "usage: " << argv[0]
			<< " <bam_path> <target> <pos> <keep_dup> <min_mapq> <min_baseq> <oxidation_error>" << endl;
		return 1;
	}

	char* path = argv[1];
	char* target = argv[2];
	int32_t pos = atol(argv[3]) - 1;
	bool keep_dup = (atoi(argv[4]) != 0);
	int min_mapq = atoi(argv[5]);
	int min_baseq = atoi(argv[6]);
	double phi = atof(argv[7]);

	cout << "global oxidation error: " << phi << endl;

	fetcher f;

	// disable quality filters
	f.qfilter.min_mapq = min_mapq;
	f.qfilter.min_baseq = min_baseq;

	if (keep_dup) {
		f.qfilter.disable_excl_flags(BAM_FDUP);
	}
	
	if (!f.open(path)) {
		cerr << "ERROR: could not open BAM file for reading: " << path << endl;
		return 1;
	}

	// fetch reads at target position
	int32_t rid = bam_name2id(f.hdr, target);
	f.fetch(rid, pos);

	size_t n = f.pile.queries.size();

	// reference base (0), alternative base (1), or other base (2)
	vector<uint8_t> bases;
	// base error probability
	vector<double> errors;
	// oxidation probability: phi if orientation is consistent, 0 otherwise
	vector<double> phis;

	vector<char> nucs;
	vector<int> quals;
	vector<char> strands;
	vector<int> read_memberships;

	bases.reserve(n);
	errors.reserve(n);
	phis.reserve(n);

	nucleotide f_ref = nuc_G, f_alt = nuc_T;
	nucleotide r_ref = nuc_C, r_alt = nuc_A;

	size_t nA = 0, nC = 0, nG = 0, nT = 0, nN = 0, nNull = 0, nDel = 0, nIns = 0;
	for (size_t i = 0; i < n; ++i) {
		bam1_t *b = f.pile.queries[i];

		nucleotide qnuc = (nucleotide) query_nucleotide(b, pos);
		if (bam_is_read1(b)) {
			if (!bam_is_rev(b)) {
				// query is on forward strand
				if (qnuc == f_ref) {
					// e.g. G>G on forward strand, read 1
					bases.push_back(0);
					phis.push_back(0);
				} else if (qnuc == f_alt) {
					// e.g. G>T on forward strand, read 1
					bases.push_back(1);
					phis.push_back(phi);
				} else if (qnuc <= nuc_T) {
					// e.g. G>A or G>C on forward strand, read 1
					// query nucleotide is a non-ref, non-alt nucleotide
					bases.push_back(2);
					phis.push_back(0);
				}
			} else {
				// query is on reverse strand: reverse complement the reference
				if (qnuc == r_ref) {
					// e.g. C>C on reverse strand, read 1
					bases.push_back(0);
					phis.push_back(0);
				} else if (qnuc == r_alt) {
					// e.g. C>A on reverse strand, read 1
					bases.push_back(1);
					phis.push_back(0);
				} else if (qnuc <= nuc_T) {
					// e.g. C>G or C>T on reverse strand, read 1
					// query nucleotide is a non-ref, non-alt nucleotide
					bases.push_back(2);
					phis.push_back(0);
				}
			}
		} else if (bam_is_read2(b)) {
			// double-check that the read is second read, in case the flag is malformed
			if (!bam_is_rev(b)) {
				// query is on forward strand
				if (qnuc == f_ref) {
					// e.g. G>G on forward strand, read 2
					bases.push_back(0);
					phis.push_back(0);
				} else if (qnuc == f_alt) {
					// e.g. G>T on forward strand, read 2
					bases.push_back(1);
					phis.push_back(0);
				} else if (qnuc <= nuc_T) {
					// query nucleotide is a non-ref, non-alt nucleotide
					// e.g. G>A or G>C on forward strand, read 2
					bases.push_back(2);
					phis.push_back(0);
				}
			} else {
				// query is on reverse strand: reverse complement the reference
				if (qnuc == r_ref) {
					// e.g. C>C on reverse strand, read 2
					bases.push_back(0);
					phis.push_back(0);
				} else if (qnuc == r_alt) {
					// e.g. C>A on reverse strand, read 2
					bases.push_back(1);
					phis.push_back(phi);
				} else if (qnuc <= nuc_T) {
					// query nucleotide is a non-ref, non-alt nucleotide
					bases.push_back(2);
					phis.push_back(0);
				}
			}
		}

		if (qnuc <= nuc_T) {
			errors.push_back(anti_phred((double)query_quality(b, pos)));
			quals.push_back(query_quality(b, pos));
			if (bam_is_rev(b)) {
				strands.push_back('-');
			} else {
				strands.push_back('+');
			}
			if (bam_is_read1(b)) {
				read_memberships.push_back(1);
			} else if (bam_is_read2(b)) {
				read_memberships.push_back(2);
			} else {
				read_memberships.push_back(0);
			}
		}

		switch (qnuc) {
		case nuc_A:
			++nA;
			nucs.push_back('A');
			break;
		case nuc_C:
			++nC;
			nucs.push_back('C');
			break;
		case nuc_G:
			++nG;
			nucs.push_back('G');
			break;
		case nuc_T:
			++nT;
			nucs.push_back('T');
			break;
		case nuc_N:
			++nN;
			break;
		case nuc_NULL:
			++nNull;
			break;
		case nuc_DEL:
			++nDel;
			break;
		case nuc_INS:
			++nIns;
			break;
		}
	}

	// print summaries
	cout << "A " << nA << endl;
	cout << "C " << nC << endl;
	cout << "G " << nG << endl;
	cout << "T " << nT << endl;
	cout << "N " << nN << endl;
	cout << "- " << nDel << endl;
	cout << "+ " << nIns << endl;
	cout << ". " << nNull << endl;
	cout << "total " << n << endl;

	cout << endl << endl;

	cout << "nuc\tbase\tqual\tstrand\tmember\terror\tphi" << endl;
	for (size_t i = 0; i < bases.size(); ++i) {
		cout 
			<< nucs[i] << '\t'
			<< (unsigned int) bases[i] << '\t'
			<< quals[i] << '\t'
			<< strands[i] << '\t'
			<< read_memberships[i] << '\t'
			<< errors[i] << '\t'
			<< phis[i] << endl;
	}

	return 0;
}
