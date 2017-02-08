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

#include "fetcher.hpp"
#include "mlat.hpp"


// region [start, end)
bool region_contains_pos(int32_t tid1, int32_t start, int32_t end, int32_t tid2, int32_t pos)
{
	if (tid1 != tid2) return false;
	if (pos >= start && pos < end) return true;
	return false;
}

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

// FIXME consider soft-clipping when extract reads?
//       bwa considers read quality and soft-clips the reads, but blat does
//       not consider read quality
//       many soft-clipped bwa alignments are not aligned by blat
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
	bool add_chr_prefix;
	int32_t n_pass_cut;

	--argc;
	if (argc == 5) {
		snv_path = argv[1];
		bam_path = argv[2];
		database = argv[3];
		out_path = argv[4];
		n_pass_cut = atoi(argv[5]);
	} else if (argc == 6) {
		snv_path = argv[1];
		bam_path = argv[2];
		database = argv[3];
		out_path = argv[4];
		n_pass_cut = atoi(argv[5]);
		db_offset = atol(argv[6]);
	} else {
		cerr << "usage: " << argv[0] <<
			" <snv.tsv> <reads.bam> <db.fasta|db.2bit> <out.vtr> <n_pass_cut> [db_offset]" << endl;
		return 1;
	}

	/// initialize database

	mlat::Database m = mlat::Database(database);

	/// initialize query

	fetcher f;

	f.qfilter.min_mapq = 5;
	f.qfilter.min_baseq = 20;
	
	f.open(bam_path);

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

	/// process each input snv

	size_t snv_i = 0;

	while (!snvf.eof()) {
		getline(snvf, line);
		if (line.empty()) break;
		if (line[0] == '#') break;

		istringstream line_stream(line);
	
		string chrom, pos_str;
		char nuc_char;

		line_stream >> chrom >> pos_str >> nuc_char;

		// set target region
		int32_t rid = bam_name2id(f.hdr, chrom.c_str());
		if (rid == -1) {
			// reference name not found: try adding chr prefix
			string chrom2 = "chr";
			chrom2 += chrom;
			rid = bam_name2id(f.hdr, chrom2.c_str());
		}

		int32_t pos = atol(pos_str.c_str()) - 1;
		nuc_t nuc = char_to_nuc(nuc_char);

		// fetch all supporting reads at query position with the mate read
		f.fetch(rid, pos, nuc, true);

		/// process each supporting read
		
		int32_t n_pass = 0;
		
		for (size_t i = 0; i < f.pile.queries.size(); ++i) {
			bam1_t *b1 = f.pile.queries[i];
			bam1_t *b2 = f.pile.mates[i];
			mlat_summary s;
			mlat_pair(m, f.pile.queries[i], f.pile.mates[i], f, rid, pos - db_offset, s);

			if (s.pass) {
				if (++n_pass > n_pass_cut) {
					// early stopping
					break;
				}
			}
		}

		if (n_pass > n_pass_cut) {
			// enough passing reads: snv passes filter
			outf << 1 << endl;
		} else {
			// snv fails filter
			outf << 0 << endl;
		}

		f.clear();

		++snv_i;
	}

	return 0;
}
