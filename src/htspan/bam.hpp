#ifndef _HTSPAN_BAM_HPP_
#define _HTSPAN_BAM_HPP_ 

/**
 * Get whether query is read 1.
 * 
 * @param b  pointer to a BAM record
 * @return true if query is read 1
 */
#define bam_is_read1(b) (((b)->core.flag & BAM_FREAD1) != 0)

/**
 * Get whether query is read 2.
 * 
 * @param b  pointer to a BAM record
 * @return true if query is read 1
 */
#define bam_is_read2(b) (((b)->core.flag & BAM_FREAD2) != 0)

#include <htslib/sam.h>
#include <string>

namespace hts {

int32_t bam_name_to_id(bam_hdr_t* hdr, const string& chrom) {
	int32_t rid = bam_name2id(hdr, chrom.c_str());
	if (rid == -1) {
		// reference name not found: try adding chr prefix
		string chrom2 = "chr";
		chrom2 += chrom;
		rid = bam_name_to_id(hdr, chrom2.c_str());
	}
	return rid;
}

}  // namespace hts

#endif  // _HTSPAN_BAM_HPP_ 

