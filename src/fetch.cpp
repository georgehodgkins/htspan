#include <iostream>
using namespace std;

#include <htslib/hts.h>
#include <htslib/bgzf.h>
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

void print_n_alignment(BGZF *fp, size_t n) {
	bam1_t *b = bam_init1();

	cout << "qname\tseq\tqual\tcigar\taux";

	for (size_t r = 0; r < n; ++r) {

		if (bam_read1(fp, b) > 0) {
			cout << bam_get_qname(b) << '\t';

			for (size_t i = 0; i < b->core.l_qseq; ++i) {
				cout << nuc_to_char(bam_seqi(bam_get_seq(b), i));
			}
			cout << '\t';

			const uint8_t* qual = bam_get_qual(b);
			for (size_t i = 0; i < b->core.l_qseq; ++i) {
				cout << char((*qual) + 33);
				++qual;
			}
			cout << '\t';

			const uint32_t* cigar = bam_get_cigar(b);
			for (size_t i = 0; i < b->core.n_cigar; ++i) {
				cout << *cigar++;
			}
			cout << '\t';

			cout << bam_get_aux(b);
		} else {
			cerr << "Error: Could not read record" << endl;
		}

	}

	bam_destroy1(b);
}

// @param ref_pos  position in the reference
void func(bam1_t *b, int ref_pos) {
	cout
		<< bam_get_qname(b) << '\t'
		<< b->core.tid << '\t' 
		<< b->core.pos << '\t' 
		<< bam_endpos(b) << '\t';

	// correct query position for preceeding cigar operations:
	// - del
	// - ref-skip
	// + soft-clip
	// + ins
	int qpos = ref_pos - b->core.pos;
	int qend = bam_endpos(b);
	// qpos will contain the nucleotide of interest in the read
	int j = 0;
	// running total of nucleotides affected by cigar operations
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

	char nuc_char = '.';

	if (j == qpos) {
		// check the next cigar operation for ins or del
		switch (bam_cigar_op(cigar[k+1])) {
				case BAM_CDEL:
					nuc_char = '-';	
					break;
				case BAM_CINS:
					nuc_char = '+';
					break;
		}
	}
	
	if (nuc_char == '.') {
		if (qpos < 0) {
			// query nucleotide is in a deletion that consumed the rest of the read
			nuc_char = '-';	
		} else if (qpos >= qend) {
			// a preceding insertion or softclip consumed the rest of the read,
			// so the read does not actually contain the sequence
			// leave the nuc_char at the default null value.
		} else {
			nuc_char = nuc_to_char(bam_seqi(bam_get_seq(b), qpos));
		}
	}

	cout << nuc_char;
	cout << endl;
}


struct plp_aux_t {
	htsFile *hf;
	bam_hdr_t *hdr;
	hts_itr_t *itr;
	hts_idx_t *idx;

	uint64_t pos0;
	int irange0;

	//bam_mate_iter_t iter;
	void *pbuffer; /* for buffered pileup */

	plp_aux_t()
		: hf(NULL), hdr(NULL), itr(NULL), idx(NULL), pbuffer(NULL) {}

	int open(const char* path) {
		clear();

		hf = hts_open(path, "rb");
		if (!hf) {
			cerr << "Error in " << __func__ << ": failed to open: " << path << endl;
			return -1;
		}
		//pos0 = bgzf_tell(hf->fp.bgzf);
		pos0 = 0;
		irange0 = 0;

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

	//int fetch(int tid, int start, int end, void *data, bam_fetch_f f) {
	int fetch_pos(int tid, int pos) {
		bam1_t *b = bam_init1();

		int start = pos, end = pos + 1;

		// seek target region
		itr = sam_itr_queryi(idx, tid, start, end);
		if (!itr) {
			cerr << "Error in " << __func__ << ": failed to seek region { tid: " 
				<< tid << ", start: " << start << ", end: " << end << " };" << endl;
		}

		cerr << "tid = " << tid << endl;
		cerr << "start = " << start << endl;
		//cerr << "seq[" << plp->qpos << "] = " << 
		//	nuc_to_char(bam_seqi(bam_get_seq(plp->b), plp->qpos)) << endl;

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

	~plp_aux_t() {
		clear();
	}
};

int fetch(const bam1_t *b, void *data) {
	return 0;
}


int plp_f(void *data, bam1_t *b) {
	cerr << "plp_f" << endl;
	plp_aux_t *p = (plp_aux_t*)data;
	int ret;
	//while ((ret = hts_itr_next(p->hf->fp.bgzf, p->itr, b, NULL)) >= 0) {
	//	func(b);
	//}
	if ((ret = hts_itr_next(p->hf->fp.bgzf, p->itr, b, NULL)) >= 0) {
		func(b, 0);
		cerr << "plp_f ret " << ret << endl;
	}
	return ret;
}

int main(int argc, char** argv) {
	char path[] = "../data/KP7-092916.bam"; 

	plp_aux_t bfile;

	cerr << "start" << endl;

	bfile.open(path);

	//cerr << h->text << endl << endl;
	cerr << "Reference sequences [" << bfile.hdr->n_targets << "]:" << endl;
	for (size_t i = 0; i < bfile.hdr->n_targets; ++i) {
		if (i == 9) {
			cerr << "..." << endl;
		} else if (i > 9 && i < bfile.hdr->n_targets - 1) {
			continue;
		}
		cerr << bfile.hdr->target_name[i] << endl;
	}
	cerr << endl;

	//print_n_alignment(bfile.hf->fp.bgzf, 3);

	// set target region
	char target[] = "chr17";
	int rid = bam_name2id(bfile.hdr, target);
	int pos = 7674420 - 1;  // 8 A, 94 G
	char seq[] = "GATGGAATCTCGCTCTGTCGCCCAGGCTGGAGTGCAGTGGCATGATCTCAGCTCACTGCAAGCTCCACCGCCCAG";
	//int pos = 7674361 - 1;   // 1 A, 2 C, 0 G, 536 T, 202 del, 15 ins
	//int pos = 7674360 - 1;    // 0 A, 759 C, 0 G, 9 T

	cerr << seq << endl;
	cerr << "chr17:7674420G>A" << endl;
	cerr << endl;
	
	bfile.fetch_pos(rid, pos);

	/*
	// seek region
	bfile.itr = sam_itr_queryi(bfile.idx, rid, start, end);

	bam_plp_t plp_itr = bam_plp_init(plp_f, &bfile);

	int max_depth = 1000;
	bam_plp_set_maxcnt(plp_itr, max_depth);

	int tid, pos, n_plp;
	const bam_pileup1_t* plp;
	while (true) {
		plp = bam_plp_auto(plp_itr, &tid, &pos, &n_plp);
		cout << "  tid = " << tid << endl;
		cout << "  pos = " << tid << endl;
		if (!plp) break;
		cout << "  qname = " << bam_get_qname(plp->b) << endl;
		cout << "  seq[" << plp->qpos << "] = " << 
			nuc_to_char(bam_seqi(bam_get_seq(plp->b), plp->qpos)) << endl;
	}
	
	bam_plp_destroy(plp_itr);
	*/
	
	return 0;
}
