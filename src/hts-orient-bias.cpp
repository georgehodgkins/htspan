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
using namespace hts;


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

	orient_bias_filter_f obfilter(nuc_G, nuc_T, phi, n);

	size_t nA = 0, nC = 0, nG = 0, nT = 0, nN = 0, nNull = 0, nDel = 0, nIns = 0;
	for (size_t i = 0; i < n; ++i) {
		bam1_t *b = f.pile.queries[i];
		obfilter.push(b, pos);

		nucleotide qnuc = (nucleotide) query_nucleotide(b, pos);
		switch (qnuc) {
		case nuc_A:
			++nA;
			break;
		case nuc_C:
			++nC;
			break;
		case nuc_G:
			++nG;
			break;
		case nuc_T:
			++nT;
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
	}  // for (size_t i = 0; i < n; ++i)

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
	for (size_t i = 0; i < obfilter.size(); ++i) {
		cout 
			<< obfilter.cnucs[i] << '\t'
			<< (int) obfilter.bases[i] << '\t'
			<< obfilter.quals[i] << '\t'
			<< obfilter.strands[i] << '\t'
			<< obfilter.members[i] << '\t'
			<< obfilter.errors[i] << '\t'
			<< obfilter.phis[i] << endl;
	}

	return 0;
}
