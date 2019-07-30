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
	
	/// vector of pointers to passing reads at the current locus
	// the pointed records are allocated and destroyed by the pileup iterator
	vector<bam1_t*> pile;

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
		plp_itr(bam_plp_init(pileup_func, this)),
		tid(0), pos(0), n(0)
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
		pile.reserve(n);
		bam_plp_set_maxcnt(plp_itr, (int) n);
	}

	size_t size() const {
		return pile.size();
	}

	/**
	 * Get reads that pass the query filter for the next site. 
	 * The pointed reads are allocated and deallocated by the 
	 * pileup iterator, so they are only valid until the next call to this
	 * function.
	 *
	 * The value of n should be checked before using the returned values;
	 * if it is 0, EOF was reached; it it is less than 0, an error occurred.
	 *
	 * @return reference to a vector of pointers to passing reads 
	 */
	const vector<bam1_t*>& next() {
		// clear buffer
		pile.clear();

		// get next pileup and update target id, position, and read depth
		const bam_pileup1_t* inner = bam_plp_auto(plp_itr, &tid, &pos, &n);

		if (!inner) {
			// if n < 0, an error occurred (not EOF)
			if (n < 0) {
				cerr << "Error in" << __func__ << ": failed to get next pileup" << endl;
			}

			return pile;
		}

		// otherwise, put each read in the returned pileup through the filter
		for (size_t i = 0; i < n; ++i) {
			if (qfilter(inner[i].b, pos)) {
				pile.push_back(inner[i].b);
			}
		}
		return pile;


	}

	/**
	 * Clear pileup.
	 */
	void clear() {
		bam_plp_reset(plp_itr);
	}

	/**
	* Return number of reads in the active pileup.
	* 0 if none loaded or -1 if last load failed.
	*/
	int curr_plp_size() const {
		return n;
	}

	/**
	* Return tid of the active pileup.
	* 0 if none loaded, unspecified if last load failed.
	*/
	int curr_plp_tid() const {
		return tid;
	}

	/**
	* Return reference position of the active pileup.
	* 0 if none loaded, unspecified if last load failed.
	*/
	int curr_plp_pos() const {
		return pos;
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
	return sam_read1(p->hf, p->hdr, b);
}

}  // namespace hts

#endif  // _HTSPAN_FETCHER_

