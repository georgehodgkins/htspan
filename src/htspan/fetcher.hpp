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

#include "frontend/simul_writer.hpp"

namespace hts {

using namespace std;

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

	/// minimum and maximum insert size for inclusion, and whether to check insert size
	int min_isize;
	int max_isize;
	bool check_isize;

	// whether to check for tandem reads, those where both the read and its mate are on the same strand
	bool excl_tandem_reads;

	// by default, exclude reads that are:
	// 1. unmapped
	// 2. secondary
	// 3. QC failed
	// 4. duplicates
	query_filter_f()
		: excl_flags(BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP), prereq_flags(0),
		min_mapq(5), min_baseq(20), min_isize(60), max_isize(600), check_isize(true), excl_tandem_reads(true)
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
		// TODO: original quality score is preferable here, use that when present
		if (b->core.qual < min_mapq) return false;

		// check that either this strand or its mate are reverse complemented, but not both
		// this check is only done if excl_tandem_reads is set to true
		bool is_rev = (b->core.flag & BAM_FREVERSE);
		bool mate_is_rev = (b->core.flag & BAM_FREVERSE);
		if (excl_tandem_reads && (is_rev == mate_is_rev)) return false;

		// check insert size, if enabled
		if (check_isize && (b->core.isize < min_isize || b->core.isize > max_isize)) return false;

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
	 * Pile is cleared before fetching.
	 *
	 * @param tid   target contig ID
	 * @param pos   target reference position
	 * @param nuc   nucleotide that the query must match at the position within
	 *              query that aligned to target reference position
	 * @param mate  whether to push mate reads too
	 * @return whether operation succeeded
	 */
	bool fetch(int32_t tid, int32_t pos, nuc_t nuc, bool mate) {
		clear();
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
	 * Close file handlers.
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

