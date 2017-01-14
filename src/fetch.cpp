#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>

enum nucleotide {
	// definitions taken from htslib
	nuc_A = 1,
	nuc_C = 2,
	nuc_G = 4,
	nuc_T = 8,
	nuc_N = 15,
	// additional definitions
	nuc_NULL = 14,
	nuc_DEL = 13,
	nuc_INS = 12
};

char nuc_to_char(uint8_t x) {
	switch (x) {
		case nuc_A: return 'A';
		case nuc_C: return 'C';
		case nuc_G: return 'G';
		case nuc_T: return 'T';
		case nuc_N: return 'N';
		case nuc_DEL: return '-';
		case nuc_INS: return '+';
		case nuc_NULL: return '.';
		default: return '.';
	}
}

uint8_t nuc_complement(uint8_t x) {
	switch (x) {
		case nuc_A: return nuc_T;
		case nuc_C: return nuc_G;
		case nuc_G: return nuc_C;
		case nuc_T: return nuc_A;
		default: return x;
	}
}

uint8_t char_to_nuc(char x) {
	switch (x) {
		case 'A': return nuc_A;
		case 'C': return nuc_C;
		case 'G': return nuc_G;
		case 'T': return nuc_T;
		case 'N': return nuc_N;
		case '-': return nuc_DEL;
		case '+': return nuc_INS;
		case '.': return nuc_NULL;
		default: return nuc_NULL;
	}
}

// check if nucleotide x and nucleotide y are equal
// if either nucleotide is N, then they are equal
bool nuc_equal(uint8_t x, uint8_t y) {
	if (x == nuc_N || y == nuc_N) return true;
	if (x == y) return true;
	return false;
}

// get the query position that aligns to a specified reference position
// if query does not align to the query position, return -1
// if query has a deletion at the specified reference position, return -2
// if query has an insertion at the specified reference position, return -3
int32_t query_position(const bam1_t *b, int32_t ref_pos) {
	// qpos will contain the nucleotide of interest in the read
	// correct query position for preceeding cigar operations:
	// - del
	// - ref-skip
	// + soft-clip
	// + ins
	int32_t qpos = ref_pos - b->core.pos;
	int32_t qend = bam_endpos(b);
	// running total of nucleotides affected by cigar operations
	int32_t j = 0, l;
	const uint32_t* cigar = bam_get_cigar(b);
	size_t n_cigar = b->core.n_cigar;
	size_t k;
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

	if (qpos >= qend) {
		// a preceding insertion or softclip consumed the rest of the read,
		// so the read does not actually contain the sequence
		// leave the nuc_c at the default null value.
		qpos = -1;
	} else if (qpos < 0) {
		// query nucleotide is in a deletion that consumed the rest of the read
		qpos = -2;	
	} else if (j == qpos) {
		// check the next cigar operation for ins or del
		switch (bam_cigar_op(cigar[k+1])) {
				case BAM_CDEL:
					qpos = -2;
					break;
				case BAM_CINS:
					qpos = -3;
					break;
		}
	}
	// otherwise, qpos points to an ordinary nucleotide within the query

	return qpos;
}

// get the query nucleotide that aligns to a specified reference position
uint8_t query_nucleotide(const bam1_t *b, int32_t ref_pos) {

	int32_t qpos = query_position(b, ref_pos);

	uint8_t nuc = nuc_NULL;
	switch (query_position(b, ref_pos)) {
		case -3:
			nuc = nuc_INS;
			break;
		case -2:
			nuc = nuc_DEL;
			break;
		case -1:
			nuc = nuc_NULL;
			break;
		default:
			if (qpos >= 0) {
				nuc = bam_seqi(bam_get_seq(b), qpos);
			} else {
				nuc = nuc_NULL;
			}
			break;
	}

	return nuc;
}

void print_query(const bam1_t *b, int32_t pos) {
	cout
		<< bam_get_qname(b) << '\t'
		<< b->core.tid << '\t' 
		<< b->core.pos << '\t' 
		<< bam_endpos(b) << '\t'
		<< nuc_to_char(query_nucleotide(b, pos))
		<< endl;
}

void print_seq(const bam1_t *b, bool original) {
	if (original && bam_is_rev(b)) {
		// query aligned to the reverse strand: the stored sequence is the reverse
		// complement of the original sequence
		// thus, reverse complement the stored sequence
		size_t n = (size_t) b->core.l_qseq;
		for (size_t i = 0; i < n; ++i) {
			cout << nuc_to_char(nuc_complement(bam_seqi(bam_get_seq(b), n - i - 1)));
		}
	} else {
		for (size_t i = 0; i < b->core.l_qseq; ++i) {
			cout << nuc_to_char(bam_seqi(bam_get_seq(b), i));
		}
	}
}

void print_seq(const bam1_t *b) {
	print_seq(b, false);
}

void print_qual(const bam1_t *b, bool original) {
	const uint8_t* qual = bam_get_qual(b);
	if (original && bam_is_rev(b)) {
		// query aligned to the reverse strand: the stored quality is reversed
		// thus, reverse the stored quality
		size_t n = (size_t) b->core.l_qseq;
		for (size_t i = 0; i < n; ++i) {
			cout << char(qual[n - i - 1] + 33);
		}
	} else {
		for (size_t i = 0; i < b->core.l_qseq; ++i) {
			cout << char((*qual) + 33);
			++qual;
		}
	}
}

void print_qual(const bam1_t *b) {
	print_qual(b, false);
}

void print_query_fasta(const bam1_t *b) {
	cout << '>' << bam_get_qname(b) << endl;
	print_seq(b, true); cout << endl;
}

void print_query_fastq(const bam1_t *b) {
	cout << '@' << bam_get_qname(b) << endl;
	print_seq(b, true); cout << endl;
	cout << '+' << endl;
	print_qual(b, true); cout << endl;
}


struct query_filter_f {
	// query will be filtered out if any of these flags is set
	int excl_flags;

	// query will be kept only if all of these flags are set
	int prereq_flags;
	
	// minimum mapping quality score for inclusion
	int min_mapq;

	// minimum base quality at query position for inclusion
	int min_baseq;

	query_filter_f()
		: excl_flags(BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP), prereq_flags(0),
		min_mapq(5), min_baseq(20)
	{}

	bool operator()(const bam1_t *b, int32_t pos) const {
		// check alignment flags
		if (prereq_flags && (b->core.flag & prereq_flags) == prereq_flags) return false;
		if (excl_flags && (b->core.flag & excl_flags) == 0) return false;

		// check mapping quality
		if (b->core.qual < min_mapq) return false;

		// check base quality at query position
		int qpos = query_position(b, pos);
		if (qpos >= 0 && bam_get_qual(b)[qpos] < min_baseq) return false;

		return true;
	}
};

struct mate_t {
	int32_t tid;
	int32_t pos;
	char *qname;

	mate_t(int32_t _tid, int32_t _pos, const char *_qname)
		: tid(_tid), pos(_pos), qname(new char[strlen(_qname)+1])
	{
		strcpy(qname, _qname);
	}

	mate_t(const mate_t& other)
		: tid(other.tid), pos(other.pos), qname(new char[strlen(other.qname)+1])
	{
		strcpy(qname, other.qname);
	}

	~mate_t() {
		delete [] qname;
	}
};

struct query_pile {

	vector<bam1_t*> queries;
	vector<bam1_t*> mates;

	// mate queries that are pending addition to `queries`
	queue<mate_t*> mate_queue;
	
	// push a query onto the pile
	// @param mate  whether to queue the mate for addition later
	bool push(const bam1_t *b, int32_t pos, uint8_t nuc, bool mate) {
		if (nuc != nuc_NULL) {
			// check if query nucleotide matches
			if (!nuc_equal(nuc, query_nucleotide(b, pos))) return false;
		}

		// copy alignment to queries
		bam1_t *b2 = bam_dup1(b);
		queries.push_back(b2);

		if (mate) {
			// add the mate of the current query to the queue
			// the mate has the same qname
			mate_queue.push( new mate_t(b2->core.mtid, b2->core.mpos, bam_get_qname(b2)) );
		}

		return true;
	}

	void clear() {
		// queries and mates should be the same size
		// but we will deallocate with separate loops to be safe
		for (size_t i = 0; i < queries.size(); ++i) {
			bam_destroy1(queries[i]);
		}
		queries.clear();
		for (size_t i = 0; i < mates.size(); ++i) {
			bam_destroy1(mates[i]);
		}
		mates.clear();
	}

	~query_pile() {
		clear();
	}
};

struct fetcher {
	// file handle
	htsFile *hf;

	// bam header
	bam_hdr_t *hdr;

	// file iterator
	hts_itr_t *itr;

	// bam index
	hts_idx_t *idx;
	
	// buffer for a single query
	bam1_t *buf;

	query_filter_f qfilter;
	query_pile pile;

	fetcher()
		: hf(NULL), hdr(NULL), itr(NULL), idx(NULL),
		buf(bam_init1())
	{
	}

	bool open(const char* path) {
		clear();

		hf = hts_open(path, "rb");
		if (!hf) {
			cerr << "Error in " << __func__ << ": failed to open: " << path << endl;
			return false;
		}

		hdr = sam_hdr_read(hf);
		if (!hdr) {
			cerr << "Error in " << __func__ << ": failed to read header for: " << path << endl;
			return false;
		}

		idx = sam_index_load(hf, path);
		if (!idx) {
			cerr << "Error in " << __func__ << ": failed to open index for: " << path << endl;
			return false;
		}

		return true;
	}

	// move iterator to target region
	bool seek(int32_t tid, int32_t start, int32_t end) {
		itr = sam_itr_queryi(idx, tid, start, end);
		if (!itr) {
			cerr << "Error in " << __func__ << ": failed to seek region { tid: " 
				<< tid << ", start: " << start << ", end: " << end << " };" << endl;
			return false;
		}
		return true;
	}

	// extract next query at current position into buffer
	bool next() {
		int ret = sam_itr_next(hf, itr, buf);
		return (ret < 0) ? false : true;
	}

	bool fetch_mates() {
		// retrieve all the mates registered in the queue from the file
		// and push them onto the mates pile in the same order,
		// inserting NULL for invalid mates
		if (pile.mate_queue.empty()) return false;

		// mates should grow to the same size as queries
		pile.mates.reserve(pile.queries.size());

		while (!pile.mate_queue.empty()) {
			bam1_t *b = NULL;
			const mate_t& m = *pile.mate_queue.front();
			if (seek(m.tid, m.pos, m.pos + 1)) {
				while (next()) {
					if (strcmp(m.qname, bam_get_qname(buf)) == 0) {
						// mate is found: make a copy
						b = bam_dup1(buf);
						// each query has only one mate
						break;
					}
				}
			}
			pile.mate_queue.pop();
			pile.mates.push_back(b);
		}

		assert(pile.queries.size() == pile.mates.size());

		return true;
	}

	bool fetch(int32_t tid, int32_t pos, uint8_t nuc, bool mate) {
		if (!seek(tid, pos, pos + 1)) return false;

		while (next()) {
			// all alignments must pass the query filter
			if (qfilter(buf, pos)) {
				// process alignment
				pile.push(buf, pos, nuc, mate);
				//print_query(buf, pos);
			}
		}

		if (mate) return fetch_mates();
		return true;
	}

	bool fetch(int32_t tid, int32_t pos) {
		return fetch(tid, pos, nuc_NULL, false);
	}

	bool fetch(const char* target_name, int32_t pos) {
		return fetch(bam_name2id(hdr, target_name), pos, nuc_NULL, false);
	}

	void clear() {
		if (idx) hts_idx_destroy(idx);
		if (itr) sam_itr_destroy(itr);
		if (hdr) bam_hdr_destroy(hdr);	
		if (hf) hts_close(hf);
	}

	~fetcher() {
		clear();
		bam_destroy1(buf);
	}
};

int main(int argc, char** argv) {
	char path[] = "../data/KP7-092916.bam"; 

	fetcher f;

	f.open(path);

	// set target region
	char target[] = "chr17";
	int32_t rid = bam_name2id(f.hdr, target);
	//int pos = 7674420 - 1;  // 8 A, 94 G
	int32_t pos = 7674361 - 1;   // 1 A, 2 C, 0 G, 536 T, 202 del, 15 ins
	//int pos = 7674360 - 1;    // 0 A, 759 C, 0 G, 9 T
	
	//f.fetch(rid, pos);
	f.fetch(rid, pos, nuc_C, true);

	for (size_t i = 0; i < f.pile.queries.size(); ++i) {
		bam1_t *b = f.pile.queries[i];
		print_query_fastq(b);
	}	
	cout << endl;
	for (size_t i = 0; i < f.pile.mates.size(); ++i) {
		bam1_t *b = f.pile.mates[i];
		print_query_fastq(b);
	}	

	return 0;
}
