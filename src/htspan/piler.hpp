#ifndef _HTSPAN_PILER_HPP_
#define _HTSPAN_PILER_HPP_

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
#include "fetcher.hpp"


namespace hts {

using namespace std;

inline int pileup_func(void* data, bam1_t* b);

/**
 * Piler of reads from BAM records.
 */
struct piler {
	/// main file handle
	htsFile *hf;

	/// pointer to BAM header
	bam_hdr_t *hdr;

	/// file iterator
	//hts_itr_t *itr;

	/// bam index
	//hts_idx_t *idx;
	
	/// buffer for a single query
	//bam1_t *buf;

	/// query filter
	query_filter_f qfilter;

	// iterator over BAM pileup
	bam_plp_t plp_itr;

	int tid;
	int pos;
	int n;

	piler()
		: hf(NULL), hdr(NULL),
		//itr(NULL), idx(NULL),
		plp_itr(bam_plp_init(pileup_func, this))
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

		/*
		idx = sam_index_load(hf, path);
		if (!idx) {
			cerr << "Error in " << __func__ << ": failed to open index for: " << path << endl;
			return false;
		}
		*/

		return true;
	}

	void reserve(int n) {
		bam_plp_set_maxcnt(plp_itr, (int) n);
	}

	size_t size() const {
		return (size_t) n;
	}

	/**
	 * Get pileup for next position.
	 *
	 * @return whether operation succeeded
	 */
	const bam_pileup1_t*  next() {
		// get next pileup and update target id, position, and read depth
		const bam_pileup1_t* inner = bam_plp_auto(plp_itr, &tid, &pos, &n);
		if (n == -1) {
			cerr << "Error in" << __func__ << ": failed to get next pileup" << endl;
			inner = NULL;
		}
		return inner;
	}

	/**
	 * Clear pileup.
	 */
	void clear() {
		bam_plp_reset(plp_itr);
	}

	/**
	 * Close file handers.
	 */
	void close() {
		//if (idx) hts_idx_destroy(idx);
		//if (itr) sam_itr_destroy(itr);
		if (hdr) {
			bam_hdr_destroy(hdr);	
			hdr = NULL;
		}
		if (hf) {
			hts_close(hf);
			hf = NULL;
		}
	}

	~piler() {
		close();

		bam_plp_destroy(plp_itr);
		plp_itr = NULL;

		//bam_destroy1(buf);
	}
};

inline int pileup_func(void* data, bam1_t* b) {
	piler* p = (piler*) data;
	int ret;
	// read next the next passing read
	do {
		ret = sam_read1(p->hf, p->hdr, b);
	} while (ret >= 0 && !p->qfilter(b, p->pos));
	return ret;	
}

}  // namespace hts

#endif  // _HTSPAN_FETCHER_

