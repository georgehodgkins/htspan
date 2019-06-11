#ifndef _HTSPAN_FAIDX_READER_HPP_
#define _HTSPAN_FAIDX_READER_HPP_

#include <iostream>
#include <vector>
#include <cassert>

#include <stdint.h>
#include <string.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "../nucleotide.hpp"


namespace hts {

using namespace std;

// TODO Do buffered reading
//      Read more sequence locally than requested and cache it

/**
 * Reader for indexed FASTA/FASTQ file
 */
struct faidx_reader {
	// main file handle
	faidx_t* fai;

	// sequence buffer
	char* seq;

	// sequence length
	int seq_len;

	faidx_reader()
	:	fai(NULL),
		seq(NULL),
		seq_len(0)
	{
	}
	
	bool open(const char* path) {
		fai = fai_load(path);
		if (fai == NULL) {
			return false;
		}
		return true;
	}

	/**
	 * Get sequence in a region
	 *
	 * Memory allocated for the sequence is managed by this class.
	 *
	 * @param tid    target id (contig index)
	 * @param start  start position (0-based)
	 * @param end    end position (0-based)
	 * @return point to sequence (NULL in case of failure)
	 */
	const char* get(int32_t tid, int32_t start, int32_t end) {
		// free previously allocated sequence
		free(seq);

		const char* c_name = faidx_iseq(fai, tid);
		seq = faidx_fetch_seq(fai, c_name, start, end, &seq_len);

		return seq;
	}

	size_t size() {
		return seq_len;
	}

	void close() {
		fai_destroy(fai);
	}

	~faidx_reader() {
		close();
		free(seq);
	}

};


}  // namespace hts

#endif  // _HTSPAN_FAIDX_READER_

