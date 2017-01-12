#include <iostream>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>

enum nucleotide {
	nuc_A = 1,
	nuc_C = 2,
	nuc_G = 4,
	nuc_T = 8,
	nuc_N = 15
};

char nuc_to_char(uint8_t x) {
	switch (x) {
		case nuc_A: return 'A';
		case nuc_C: return 'C';
		case nuc_G: return 'G';
		case nuc_T: return 'T';
		default: return 'N';
	}
}

uint8_t char_to_nuc(char x) {
	switch (x) {
		case 'A': return nuc_A;
		case 'C': return nuc_C;
		case 'G': return nuc_G;
		case 'T': return nuc_T;
		default: return nuc_N;
	}
}

// @param ref_pos  position in the reference
void func(bam1_t *b, int ref_pos) {
	cout
		<< bam_get_qname(b) << '\t'
		<< b->core.tid << '\t' 
		<< b->core.pos << '\t' 
		<< bam_endpos(b) << '\t';

	// qpos will contain the nucleotide of interest in the read
	// correct query position for preceeding cigar operations:
	// - del
	// - ref-skip
	// + soft-clip
	// + ins
	int qpos = ref_pos - b->core.pos;
	int qend = bam_endpos(b);
	// running total of nucleotides affected by cigar operations
	int j = 0;
	int n_cigar = b->core.n_cigar;
	const uint32_t* cigar = bam_get_cigar(b);
	int l, k;
	for (k = 0; k < n_cigar; ++k) {
		l = bam_cigar_oplen(cigar[k]);
		j += l;
		switch (bam_cigar_op(cigar[k])) {
				case BAM_CMATCH:
				case BAM_CEQUAL:
				case BAM_CDIFF:
					break;
				case BAM_CREF_SKIP:
				case BAM_CDEL:
					qpos -= l;
					break;
				case BAM_CSOFT_CLIP:
				case BAM_CINS:
					qpos += l;
					break;
				// BAM_CPAD ????
		}
		if (j >= qpos) break;
	}

	char nuc_c = '.';

	if (j == qpos) {
		// check the next cigar operation for ins or del
		switch (bam_cigar_op(cigar[k+1])) {
				case BAM_CDEL:
					nuc_c = '-';	
					break;
				case BAM_CINS:
					nuc_c = '+';
					break;
		}
	}
	
	if (nuc_c == '.') {
		if (qpos < 0) {
			// query nucleotide is in a deletion that consumed the rest of the read
			nuc_c = '-';	
		} else if (qpos >= qend) {
			// a preceding insertion or softclip consumed the rest of the read,
			// so the read does not actually contain the sequence
			// leave the nuc_c at the default null value.
		} else {
			nuc_c = nuc_to_char(bam_seqi(bam_get_seq(b), qpos));
		}
	}

	cout << nuc_c;
	cout << endl;
}

struct bamfile {
	htsFile *hf;
	bam_hdr_t *hdr;
	hts_itr_t *itr;
	hts_idx_t *idx;

	bamfile()
		: hf(NULL), hdr(NULL), itr(NULL), idx(NULL) {}

	int open(const char* path) {
		clear();

		hf = hts_open(path, "rb");
		if (!hf) {
			cerr << "Error in " << __func__ << ": failed to open: " << path << endl;
			return -1;
		}

		hdr = sam_hdr_read(hf);
		if (!hdr) {
			cerr << "Error in " << __func__ << ": failed to read header for: " << path << endl;
			return -1;
		}

		idx = sam_index_load(hf, path);
		if (!idx) {
			cerr << "Error in " << __func__ << ": failed to open index for: " << path << endl;
			return -1;
		}

		return 0;
	}

	int fetch_pos(int tid, int pos) {
		bam1_t *b = bam_init1();

		int start = pos, end = pos + 1;

		// seek target region
		itr = sam_itr_queryi(idx, tid, start, end);
		if (!itr) {
			cerr << "Error in " << __func__ << ": failed to seek region { tid: " 
				<< tid << ", start: " << start << ", end: " << end << " };" << endl;
		}

		// extract alignments to b at itr
		while (sam_itr_next(hf, itr, b) >= 0) {
			func(b, pos);
		}

		bam_destroy1(b);

		return 0;
	}

	int fetch_pos(const char* target_name, int pos) {
		return fetch_pos(bam_name2id(hdr, target_name), pos);
	}

	void clear() {
		if (idx) hts_idx_destroy(idx);
		if (itr) sam_itr_destroy(itr);
		if (hdr) bam_hdr_destroy(hdr);	
		if (hf) hts_close(hf);
	}

	~bamfile() {
		clear();
	}
};

int main(int argc, char** argv) {
	char path[] = "../data/KP7-092916.bam"; 

	bamfile bfile;

	bfile.open(path);

	// set target region
	char target[] = "chr17";
	int rid = bam_name2id(bfile.hdr, target);
	int pos = 7674420 - 1;  // 8 A, 94 G
	char seq[] = "GATGGAATCTCGCTCTGTCGCCCAGGCTGGAGTGCAGTGGCATGATCTCAGCTCACTGCAAGCTCCACCGCCCAG";
	//int pos = 7674361 - 1;   // 1 A, 2 C, 0 G, 536 T, 202 del, 15 ins
	//int pos = 7674360 - 1;    // 0 A, 759 C, 0 G, 9 T
	
	bfile.fetch_pos(rid, pos);
	
	return 0;
}
