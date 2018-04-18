#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
#include <vector>
#include <cassert>
#include <cstdlib>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "htspan/fetcher.hpp"
using namespace hts;


void read_snv_record(const string& line, const fetcher& f, int32_t& rid, 
		int32_t& pos, nuc_t& ref, nuc_t& alt) {
	
	istringstream stream(line);

	string chrom, pos_str;
	char ref_c, alt_c;
	stream >> chrom >> pos_str >> ref_c >> alt_c;

	rid = bam_name2id(f.hdr, chrom.c_str());
	if (rid == -1) {
		// reference name not found: try adding chr prefix
		string chrom2 = "chr";
		chrom2 += chrom;
		rid = bam_name2id(f.hdr, chrom2.c_str());
	}

	pos = atol(pos_str.c_str()) - 1;
	ref = char_to_nuc(ref_c);
	alt = char_to_nuc(alt_c);
}

template<class C, class T>
int process(basic_istream<C, T>& snvf, basic_ostream<C, T>& outf, fetcher& f, bool echo_snv) {

	/// initialize input file

	string line;
	// discard header line
	getline(snvf, line);

	/// initialize output file

	// write header line
	if (echo_snv) {
		outf << "chrom\tpos\tref\talt\t";
	}
	outf << "ref_count\talt_count" << endl;

	/// process each input snv

	size_t snv_i = 0;

	while (!snvf.eof()) {
		getline(snvf, line);
		if (line.empty()) break;
		if (line[0] == '#') continue;

		int32_t rid, pos;
		nuc_t ref, alt;

		read_snv_record(line, f, rid, pos, ref, alt);

		// fetch all reads at query position without mate
		f.fetch(rid, pos);

		// tally the alleles at the query position
		size_t n = f.pile.queries.size();
		size_t nA = 0, nC = 0, nG = 0, nT = 0, nN = 0, nNull = 0, nDel = 0, nIns = 0;
		for (size_t i = 0; i < n; ++i) {
			bam1_t *b = f.pile.queries[i];
			switch (query_nucleotide(b, pos)) {
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
		}

		int32_t ref_count;
		switch (ref) {
		case nuc_A:
			ref_count = nA;
			break;
		case nuc_C:
			ref_count = nC;
			break;
		case nuc_G:
			ref_count = nG;
			break;
		case nuc_T:
			ref_count = nT;
			break;
		default:
			ref_count = 0;
			break;
		}

		int32_t alt_count;
		switch (alt) {
		case nuc_A:
			alt_count = nA;
			break;
		case nuc_C:
			alt_count = nC;
			break;
		case nuc_G:
			alt_count = nG;
			break;
		case nuc_T:
			alt_count = nT;
			break;
		default:
			alt_count = 0;
			break;
		}

		// write output line
		if (echo_snv) {
			outf << line << '\t';
		}
		outf << ref_count << '\t' << alt_count << endl;

		f.clear();

		++snv_i;
	}

	return 0;
}


int main(int argc, char** argv) {
	--argc;

	if (argc != 7) {
		cerr << "usage: " << argv[0]
			<< " <snv4.tsv> <reads.bam> <cov.tsv> <keep_dup> <min_mapq> <min_baseq> <echo_snv>" << endl << endl;
		cerr << "Format specification of the SNV4 file: " << endl;
		cerr << "snv4file ::=  header [TAB record]+" << endl;
		cerr << "header   ::=  'chrom' TAB 'pos' TAB 'ref' TAB 'alt' CRLF" << endl;
		cerr << "record   ::=  chrom TAB pos TAB ref TAB alt CRLF" << endl;
		cerr << "chrom    ::=  ('chr' | '') ([0-9]+ | 'X' | 'Y' | 'MT')" << endl;
		cerr << "pos      ::=  [0-9]+" << endl;
		cerr << "alt      ::=  nucl" << endl;
		cerr << "ref      ::=  nucl" << endl;
		cerr << "nucl     ::=  'A' | 'C' | 'G' | 'T'" << endl;
		return 1;
	}

	char* snv_path = argv[1];
	char* bam_path = argv[2];
	char* out_path = argv[3];

	bool keep_dup = (atoi(argv[4]) != 0);
	int min_mapq = atoi(argv[5]);
	int min_baseq = atoi(argv[6]);
	bool echo_snv = (atoi(argv[7]) != 0);


	/// initialize query

	fetcher f;

	f.qfilter.min_mapq = min_mapq;
	f.qfilter.min_baseq = min_baseq;

	if (keep_dup) {
		f.qfilter.disable_excl_flags(BAM_FDUP);
	}
	
	if (!f.open(bam_path)) {
		cerr << "ERROR: could not open BAM file for reading: " << bam_path << endl;
		return 1;
	}

	/// main process
	
	istream *pfin;
	ostream *pfout;

	// use stdin or stdout if file name is "-"

	ifstream fin;
	if (strcmp(snv_path, "-") != 0) {
		fin.open(snv_path);
		if (!fin.is_open()) {
			cerr << "ERROR: could not open input file" << endl;
			return 1;
		}
		pfin = &fin;
	} else {
		pfin = &cin;
	}

	ofstream fout;
	if (strcmp(out_path, "-") != 0) {
		fout.open(out_path);
		if (!fout.is_open()) {
			cerr << "ERROR: could not open output file" << endl;
			return 1;
		}
		pfout = &fout;
	} else {
		pfout = &cout;
	}
	
	int ret = process(*pfin, *pfout, f, echo_snv);

	fin.close();
	fout.close();

	return ret;
}
