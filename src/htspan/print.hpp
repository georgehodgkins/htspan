#ifndef _HTSPAN_PRINT_HPP_
#define _HTSPAN_PRINT_HPP_

#include <iostream>

#include "fetcher.hpp"

namespace hts {

using namespace std;

void print_seq(const bam1_t *b, bool original) {
	if (b == NULL) return;
	size_t n = (size_t) b->core.l_qseq;
	if (original && bam_is_rev(b)) {
		// query aligned to the reverse strand: the stored sequence is the reverse
		// complement of the original sequence
		// thus, reverse complement the stored sequence
		for (size_t i = 0; i < n; ++i) {
			cout << nuc_to_char(nuc_complement(bam_seqi(bam_get_seq(b), n - i - 1)));
		}
	} else {
		for (size_t i = 0; i < n; ++i) {
			cout << nuc_to_char(bam_seqi(bam_get_seq(b), i));
		}
	}
}

void print_seq(const bam1_t *b) {
	print_seq(b, false);
}

void print_qual(const bam1_t *b, bool original) {
	if (b == NULL) return;
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
	if (b == NULL) return;
	cout << '>' << bam_get_qname(b) << endl;
	print_seq(b, true); cout << endl;
}

void print_query_fastq(const bam1_t *b) {
	if (b == NULL) return;
	cout << '@' << bam_get_qname(b) << endl;
	print_seq(b, true); cout << endl;
	cout << '+' << endl;
	print_qual(b, true); cout << endl;
}

void print_query(const bam1_t *b, int32_t pos) {
	if (b == NULL) return;
	cout
		<< bam_get_qname(b) << '\t'
		<< b->core.tid << '\t' 
		<< b->core.pos << '\t' 
		<< bam_endpos(b) << '\t'
		<< nuc_to_char(query_nucleotide(b, pos))
		<< endl;
}

void print_references(fetcher& f, size_t n) {
	cout << "reference sequences [" << f.hdr->n_targets << "]:" << endl;
	const int total = f.hdr->n_targets;
	for (size_t i = 0; i < total; ++i) {
		if (i == n - 2 && n < total) {
			cout << "..." << endl;
		} else if (i > n - 2 && i < f.hdr->n_targets - 1) {
			// print nothing
		} else {
			cout << f.hdr->target_name[i] << endl;
		}
	}
	cout << endl;
}

void print_n_alignment(fetcher& f, size_t n) {
	bam1_t *b = bam_init1();

	cout << "qname\tseq\tqual\tcigar\taux" << endl;

	for (size_t r = 0; r < n; ++r) {

		if (bam_read1(f.hf->fp.bgzf, b) > 0) {
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

			cout << endl;
		} else {
			cerr << "Error: Could not read record" << endl;
		}

	}
	cout << endl;

	bam_destroy1(b);
}


}  // namespace hts


#endif  // _HTSPAN_PRINT_HPP_

