#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <string>
#include <algorithm>
#include <limits>
#include <fstream>
#include <sstream>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <mlat.hpp>

#include "htspan/fetcher.hpp"
using namespace hts;


// region [start, end)
bool region_contains_pos(int32_t tid1, int32_t start, int32_t end, int32_t tid2, int32_t pos)
{
	if (tid1 != tid2) return false;
	if (pos >= start && pos < end) return true;
	return false;
}

/*
template <typename T>
T min4(T x1, T x2, T x3, T x4) {
	return std::min(std::min(std::min(x1, x2), x3), x4);
}

template <typename T>
T max4(T x1, T x2, T x3, T x4) {
	return std::max(std::max(std::max(x1, x2), x3), x4);
}
*/

struct mlat_summary
{
	// mlat alignment summary for read1 and read2
	//
	// index of best alignments in mlat results
	size_t idx1, idx2;

	// target coordinates of the best alignment
	int32_t tid1, start1, end1, tid2, start2, end2;

	// differences in mlat scores
	int sdiff1, sdiff2;

	// fragment size
	int32_t frag_size;

	bool pass;

	mlat_summary()
		: tid1(0), start1(0), end1(0), sdiff1(0),
		tid2(0), start2(0), end2(0), sdiff2(0),
		idx1(0), idx2(0), frag_size(0), 
		pass(true) {}
};

void write_mlat_summary_header(ofstream& fout, bool show_coords) {
	if (show_coords) {
		fout <<
			"tid1\t"
			"start1\t"
			"end1\t"
			"tid2\t"
			"start2\t"
			"end2\t";
	}
	fout <<
		"idx1\t"
		"idx2\t"
		"sdiff1\t"
		"sdiff2\t"
		"frag_size\t"
		"pass";
}

void write_mlat_summary(ofstream& fout, const mlat_summary& s, int32_t offset, bool show_coords) {
	if (show_coords) {
		fout
			<< s.tid1+1 << '\t'
			<< s.start1+1 + offset << '\t'
			<< s.end1 + offset << '\t'
			<< s.tid2+1 << '\t'
			<< s.start2+1 + offset << '\t'
			<< s.end2 + offset << '\t';
	}
	fout
		<< s.idx1+1 << '\t'
		<< s.idx2+1 << '\t'
		<< s.sdiff1 << '\t'
		<< s.sdiff2 << '\t'
		<< s.frag_size << '\t'
		<< s.pass;
}

// FIXME consider soft-clipping when extract reads?
// b2 may be NULL
void mlat_pair(mlat::Database& m, bam1_t *b1, bam1_t *b2, const fetcher& f, int32_t tid, int32_t pos, mlat_summary& out) {

	// re-align the first read of the pair
	char *seq1 = bam_seq_cstr(b1, true);
	mlat::Result res1(m.search(seq1));

	// find best alignment that contains query position
	// (results are sorted by desending score)
	size_t valid_i = 0;
	while (valid_i < res1.size()) {
		out.tid1 = bam_name2id(f.hdr, res1[valid_i].tName);
		out.start1 = res1[valid_i].tStart;
		out.end1 = res1[valid_i].tEnd;
		if (region_contains_pos(out.tid1, out.start1, out.end1, tid, pos)) {
			break;
		}
		++valid_i;
	}
	out.idx1 = valid_i;
	// index of the alignment with the next best score
	size_t next_i = (valid_i == 0) ? 1 : 0;

	if (valid_i < res1.size()) {

		out.sdiff1 = res1[valid_i].score;
		if (next_i < res1.size()) {
			out.sdiff1 -= res1[next_i].score;
		}

		if (b2 != NULL) {
			//assert( strcmp(bam_get_qname(b1), bam_get_qname(b2)) == 0 );

			// re-align the second read of the pair
			char *seq2 = bam_seq_cstr(b2, true);
			mlat::Result res2(m.search(seq2));
			delete [] seq2;

			// find alignment of the second read closest to the alignment of the first read
			size_t argmin_j = 0;
			int32_t min_frag_size = std::numeric_limits<int32_t>::max();
			for (size_t j = 0; j < res2.size(); ++j) {
				out.tid2 = bam_name2id(f.hdr, res2[j].tName);
				out.start2 = res2[j].tStart;
				out.end2 = res2[j].tEnd;

				if (out.tid2 == out.tid1) {
					// since BLAT only searches the + strand of the target (it searches the query and the
					// reverse complement of the query), tStart is always the leftmost
					// coordinate and tEnd is always the rightmost coordinate
					// therefore, we only need to compare tStart and tEnd to find the
					// start and end of the read-pair fragment
					int32_t frag_start = std::min(out.start1, out.start2);
					int32_t frag_end = std::max(out.end1, out.end2);
					int32_t frag_size = frag_end - frag_start;
					if (frag_size < min_frag_size) {
						min_frag_size = frag_size;
						argmin_j = j;
					}
				}
			}

			if (min_frag_size != std::numeric_limits<int32_t>::max()) {
				out.idx2 = argmin_j;
				out.frag_size = min_frag_size;
				size_t next_j = (argmin_j == 0) ? 1 : 0;
				out.sdiff2 = res2[argmin_j].score;
				if (next_j < res2.size()) {
					out.sdiff2 -= res2[next_j].score;
				}
			}
		}

	} else {
		// none of the alignment of the first read contains the position
		out.pass = false;
	}

	if (out.sdiff1 > 0 || ((out.sdiff1 == 0) && out.sdiff2 > 0)) {
		out.pass = true;
	} else {
		out.pass = false;
		cerr << "INFO: read pair failed: " << bam_get_qname(b1) << ' ' << seq1 << endl;
	}

	delete [] seq1;
}

int main(int argc, char** argv) {

	char *snv_path, *bam_path, *database, *out_path;
	int32_t db_offset;
	bool add_chr_prefix, show_coords;

	--argc;
	if (argc == 4) {
		snv_path = argv[1];
		bam_path = argv[2];
		database = argv[3];
		out_path = argv[4];
		db_offset = 0;
		add_chr_prefix = 0;
		show_coords = 0;
	} else if (argc == 7) {
		snv_path = argv[1];
		bam_path = argv[2];
		database = argv[3];
		out_path = argv[4];
		db_offset = atol(argv[5]);
		add_chr_prefix = (atoi(argv[6]) != 0);
		show_coords = (atoi(argv[7]) != 0);
	} else {
		cerr << "usage: " << argv[0] <<
			" <snv3.tsv> <reads.bam> <db.fasta|db.2bit> <out.tsv> <db_offset> <add_chr_prefix> <show_coords>" << endl << endl;
		cerr << "Format specification of the SNV3 file: " << endl;
		cerr << "snv3file ::=  header [TAB record]+" << endl;
		cerr << "header   ::=  'chrom' TAB 'pos' TAB 'alt' CRLF" << endl;
		cerr << "record   ::=  chrom TAB pos TAB alt CRLF" << endl;
		cerr << "chrom    ::=  ('chr' | '') ([0-9]+ | 'X' | 'Y' | 'MT')" << endl;
		cerr << "pos      ::=  [0-9]+" << endl;
		cerr << "alt      ::=  'A' | 'C' | 'G' | 'T'" << endl;
		return 1;
	}

	/// initialize database

	mlat::Database m = mlat::Database(database);

	/// initialize query

	fetcher f;

	f.qfilter.min_mapq = 5;
	f.qfilter.min_baseq = 20;
	
	if (!f.open(bam_path)) {
		cerr << "ERROR: could not open BAM file for reading: " << bam_path << endl;
		return 1;
	}

	/// initialize input file

	ifstream snvf(snv_path);

	if (!snvf.is_open()) {
		cerr << "ERROR: could not open for read: " << snv_path << endl;
		return 1;
	}

	string line;
	// discard header line
	getline(snvf, line);

	/// initialize output file

	ofstream outf(out_path);

	if (!outf.is_open()) {
		cerr << "ERROR: could not open for write: " << out_path << endl;
		return 1;
	}

	outf << "query\treadpair\t";
	write_mlat_summary_header(outf, show_coords);
	outf << endl;

	/// process each input snv

	size_t snv_i = 0;

	while (!snvf.eof()) {
		getline(snvf, line);
		if (line.empty()) break;
		if (line[0] == '#') break;

		istringstream line_stream(line);
	
		string chrom, target = "", pos_str;
		char nuc_char;

		line_stream >> chrom >> pos_str >> nuc_char;

		if (add_chr_prefix) {
			if (chrom.substr(0, 3) != "chr") {
				target += "chr";
				target += chrom;
			} else {
				target = chrom;
			}
		} else {
			target = chrom;
		}

		// set target region
		int32_t rid = bam_name2id(f.hdr, target.c_str());
		int32_t pos = atol(pos_str.c_str()) - 1;
		nuc_t nuc = char_to_nuc(nuc_char);

		// fetch all supporting reads at query position with the mate read
		f.fetch(rid, pos, nuc, true);

		/// process each supporting read
		
		for (size_t i = 0; i < f.pile.queries.size(); ++i) {
			bam1_t *b1 = f.pile.queries[i];
			bam1_t *b2 = f.pile.mates[i];
			mlat_summary s;
			mlat_pair(m, f.pile.queries[i], f.pile.mates[i], f, rid, pos - db_offset, s);

			outf << snv_i+1 << '\t' << i+1 << '\t';
			write_mlat_summary(outf, s, db_offset, show_coords);
			outf << endl;
		}

		f.clear();

		++snv_i;
	}

	return 0;
}
