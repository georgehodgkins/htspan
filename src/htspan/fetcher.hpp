#ifndef _HTSPAN_FETCHER_HPP_
#define _HTSPAN_FETCHER_HPP_

#include <iostream>
#include <queue>
#include <vector>
#include <cassert>

#include <stdint.h>
#include <string.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "bam.hpp"
#include "nucleotide.hpp"


namespace hts {

using namespace std;

/**
 * Get the query position that aligns to a specified reference position.
 *
 * Query position is a position within the query read.
 * If query does not align to the query position, return -1.
 * If query has a deletion at the specified reference position, return -2.
 * If query has an insertion at the specified reference position, return -3.
 * 
 * @param b        pointer to BAM record
 * @param ref_pos  reference position
 * @return query position
 */
int32_t query_position(const bam1_t *b, int32_t ref_pos) {
	// qpos will contain the nucleotide of interest in the read
	// correct query position for preceeding cigar operations:
	// - del
	// - ref-skip
	// + soft-clip
	// + ins
	// cigar string is in the same orientation as the stored reads
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

/**
 * Get the query nucleotide that aligns to a specified reference position.
 *
 * The nucleotide is reported with respect to the reference forward 
 * strand, so reads aligned to the reverse strand will have been 
 * reverse-complemented,
 *
 * @param b        pointer to the BAM record
 * @param ref_pos  reference position
 * @return nucleotide enum
 */
nuc_t query_nucleotide(const bam1_t *b, int32_t ref_pos) {
	int32_t qpos = query_position(b, ref_pos);

	uint8_t nuc = nuc_NULL;
	switch (qpos) {
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

/**
 * Get the query quality that aligns to a specified reference position.
 *
 * @param b        pointer to the BAM record
 * @param ref_pos  reference position
 * @return Phred base quality score
 */
uint8_t query_quality(const bam1_t *b, int32_t ref_pos) {
	int32_t qpos = query_position(b, ref_pos);
	if (qpos >= 0) {
		return bam_get_qual(b)[qpos];
	}
	return 0;
}

/**
 * Get the sequence from a BAM record.
 *
 * @param b         pointer to BAM record
 * @param s         output string
 * @param original  whether to return the original read sequence
 */
void bam_seq_str(const bam1_t *b, string& s, bool original) {
	if (b == NULL) return;
	size_t n = (size_t) b->core.l_qseq;
	s.resize(n);
	if (original && bam_is_rev(b)) {
		// query aligned to the reverse strand: the stored sequence is the reverse
		// complement of the original sequence
		// thus, reverse complement the stored sequence
		for (size_t i = 0; i < n; ++i) {
			s[i] = nuc_to_char(nuc_complement(bam_seqi(bam_get_seq(b), n - i - 1)));
		}
	} else {
		for (size_t i = 0; i < n; ++i) {
			s[i] = nuc_to_char(bam_seqi(bam_get_seq(b), i));
		}
	}
}

void bam_seq_str(const bam1_t *b, string& s) {
	bam_seq_str(b, s, false);
}

/**
 * Get the sequence from a BAM record.
 *
 * @param b         pointer to BAM record
 * @param original  whether to return the original read sequence
 * @return pointer to newly allocated C-string (ownership is moved to caller)
 */
char *bam_seq_cstr(const bam1_t *b, bool original) {
	if (b == NULL) return NULL;
	size_t n = (size_t) b->core.l_qseq;
	char* s = new char[n + 1];
	if (original && bam_is_rev(b)) {
		// query aligned to the reverse strand: the stored sequence is the reverse
		// complement of the original sequence
		// thus, reverse complement the stored sequence
		for (size_t i = 0; i < n; ++i) {
			s[i] = nuc_to_char(nuc_complement(bam_seqi(bam_get_seq(b), n - i - 1)));
		}
	} else {
		for (size_t i = 0; i < n; ++i) {
			s[i] = nuc_to_char(bam_seqi(bam_get_seq(b), i));
		}
	}
	s[n] = '\0';
	return s;
}

char *bam_seq_cstr(const bam1_t *b) {
	return bam_seq_cstr(b, false);
}

/**
 * Functor for query filtering.
 */
struct query_filter_f {
	/// query will be filtered out if any of these flags is set
	int excl_flags;

	/// query will be kept only if all of these flags are set
	int prereq_flags;
	
	/// minimum mapping quality score for inclusion
	int min_mapq;

	/// minimum base quality at query position for inclusion
	int min_baseq;

	// by default, exclude reads that are:
	// 1. unmapped
	// 2. secondary
	// 3. QC failed
	// 4. duplicates
	query_filter_f()
		: excl_flags(BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP), prereq_flags(0),
		min_mapq(5), min_baseq(20)
	{}

	void enable_excl_flags(int flags) {
		excl_flags |= flags;
	}

	void disable_excl_flags(int flags) {
		excl_flags &= ~flags;
	}

	void enable_prereq_flags(int flags) {
		prereq_flags |= flags;
	}

	void disable_prereq_flags(int flags) {
		prereq_flags &= ~flags;
	}

	bool operator()(const bam1_t *b, int32_t pos) const {
		// check alignment flags
		if (prereq_flags && (b->core.flag & prereq_flags) != prereq_flags) return false;
		if (excl_flags && (b->core.flag & excl_flags) != 0) return false;

		// check mapping quality
		if (b->core.qual < min_mapq) return false;

		// check base quality at query position
		int qpos = query_position(b, pos);
		if (qpos >= 0 && bam_get_qual(b)[qpos] < min_baseq) return false;

		return true;
	}
};

/**
 * Mate query.
 */
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

/**
 * Piles of query reads.
 *
 * Stores the BAM records of paired-end reads.
 */
struct query_pile {

	/// BAM record of first read
	vector<bam1_t*> queries;

	/// BAM record of mate/second read
	vector<bam1_t*> mates;

	/// mate queries that are pending addition to `queries`
	queue<mate_t*> mate_queue;
	
	/** 
	 * Push a query onto the pile.
	 *
	 * Query read is only pushed if the query nucleotide that aligned to
	 * the target reference position matches the specified nucleotide.
	 *
	 * @param b     pointer to BAM record
	 * @param pos   reference position
	 * @param nuc   nucleotide that the query must match at the position within
	 *              query that aligned to target reference position
	 * @param mate  whether to queue the mate for addition later
	 * @return whether operation succeeded
	 */
	bool push(const bam1_t *b, int32_t pos, nuc_t nuc, bool mate) {
		if (nuc != nuc_NULL) {
			// check if query nucleotide matches
			if (!nuc_equal(nuc, query_nucleotide(b, pos))) return false;
		}

		// copy alignment to queries
		bam1_t *b2 = bam_dup1(b);
		queries.push_back(b2);

		if (mate) {
			// Add the mate of the current query to the queue
			// the mate has the same qname.
			// It is more efficient to access reads that are close together and
			// avoid IO seeks.
			// Caller is repsonsible for actually pushing mate onto mates pile.
			mate_queue.push( new mate_t(b2->core.mtid, b2->core.mpos, bam_get_qname(b2)) );
		}

		return true;
	}

	/**
	 * Dellocate all data.
	 */
	void clear() {
		// Queries and mates may not be same size, depending on whether caller
		// pushed reads onto reads; therefore we deallocate using separate loops
		for (size_t i = 0; i < queries.size(); ++i) {
			bam_destroy1(queries[i]);
		}
		queries.clear();
		for (size_t i = 0; i < mates.size(); ++i) {
			bam_destroy1(mates[i]);
		}
		mates.clear();

		while (!mate_queue.empty()) {
			delete mate_queue.front();
			mate_queue.pop();
		}
	}

	~query_pile() {
		clear();
	}
};

/**
 * Fetcher of reads from BAM records.
 */
struct fetcher {
	/// main file handle
	htsFile *hf;

	/// pointer to BAM header
	bam_hdr_t *hdr;

	/// file iterator
	hts_itr_t *itr;

	/// bam index
	hts_idx_t *idx;
	
	/// buffer for a single query
	bam1_t *buf;

	/// query filter
	query_filter_f qfilter;

	/// pile of query reads
	query_pile pile;

	fetcher()
		: hf(NULL), hdr(NULL), itr(NULL), idx(NULL),
		buf(bam_init1())
	{
	}

	/**
	 * Open BAM file and its index.
	 *
	 * @param path  file path to BAM file
	 * @return whether operation succeeded
	 */
	bool open(const char* path) {
		close();

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

	/**
	 * Move iterator to target region.
	 *
	 * @param tid    target contig ID (e.g. chromosome name of reference)
	 * @param start  start position of region
	 * @param end    end position of region
	 * @return whether operation succeeded
	 */
	bool seek(int32_t tid, int32_t start, int32_t end) {
		itr = sam_itr_queryi(idx, tid, start, end);
		if (!itr) {
			cerr << "Error in " << __func__ << ": failed to seek region { tid: " 
				<< tid << ", start: " << start << ", end: " << end << " };" << endl;
			return false;
		}
		return true;
	}

	/**
	 * Extract next query at current position into buffer.
	 *
	 * @return whether operation succeeded
	 */
	bool next() {
		int ret = sam_itr_next(hf, itr, buf);
		return (ret < 0) ? false : true;
	}

	/**
	 * Fetch the mate reads and push onto pile.
	 *
	 * @return whether operation succeeded
	 */
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

	/**
	 * Fetch read-pairs at position and push onto pile.
	 *
	 * @param tid   target contig ID
	 * @param pos   target reference position
	 * @param nuc   nucleotide that the query must match at the position within
	 *              query that aligned to target reference position
	 * @param mate  whether to push mate reads too
	 * @return whether operation succeeded
	 */
	bool fetch(int32_t tid, int32_t pos, nuc_t nuc, bool mate) {
		if (!seek(tid, pos, pos + 1)) return false;

		while (next()) {
			// all alignments must pass the query filter
			if (qfilter(buf, pos)) {
				// process alignment
				pile.push(buf, pos, nuc, mate);
			}
		}

		if (mate) return fetch_mates();
		return true;
	}

	/**
	 * Fetch reads at target position.
	 *
	 * Mate reads are not fetched.
	 *
	 * @param tid  target contig ID
	 * @param pos  target reference position
	 * @return whether operation succeeded
	 */
	bool fetch(int32_t tid, int32_t pos) {
		return fetch(tid, pos, nuc_NULL, false);
	}

	/**
	 * Fetch reads at target position.
	 *
	 * Mate reads are not fetched.
	 *
	 * @param target_name  target name (e.g. chromosome name of reference)
	 * @param pos          target reference position
	 * @return whether operation succeeded
	 */
	bool fetch(const char* target_name, int32_t pos) {
		return fetch(bam_name2id(hdr, target_name), pos, nuc_NULL, false);
	}

	/**
	 * Clear read pile.
	 */
	void clear() {
		pile.clear();
	}

	/**
	 * Close file handers.
	 */
	void close() {
		if (idx) hts_idx_destroy(idx);
		if (itr) sam_itr_destroy(itr);
		if (hdr) bam_hdr_destroy(hdr);	
		if (hf) hts_close(hf);
	}

	~fetcher() {
		close();
		bam_destroy1(buf);
	}
};

}  // namespace hts

#endif  // _HTSPAN_FETCHER_

