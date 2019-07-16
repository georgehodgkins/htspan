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

#include "nucleotide.hpp"

namespace hts {

using namespace std;

int32_t bam_name_to_id(bam_hdr_t* hdr, const string& chrom, bool first_pass = true) {
	int32_t rid = bam_name2id(hdr, chrom.c_str());
	if (rid == -1 && first_pass) {
		// reference name not found: try adding chr prefix
		string chrom2 = "chr";
		chrom2 += chrom;
		rid = bam_name_to_id(hdr, chrom2.c_str(), false);
	}
	return rid;
}

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

}  // namespace hts

#endif  // _HTSPAN_BAM_HPP_ 

