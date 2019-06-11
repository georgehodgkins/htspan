#ifndef _HTSPAN_PRINT_HPP_
#define _HTSPAN_PRINT_HPP_
//TODO move to src/bin

#include <iostream>

#include "fetcher.hpp"

namespace hts {

using namespace std;

/**
* Print the sequence stored in a given BAM record.
* 
* @param b Pointer to BAM record containing sequence to be printed.
* @param original Whether the original (non reverse-complemented)
*	seq should be returned for a reverse strand read. [false]
* @return none (prints the sequence in the record to stdout)
*/
void print_seq(const bam1_t *b, bool original = false) {
	if (b == NULL) return;
	int n = b->core.l_qseq;
	if (original && bam_is_rev(b)) {
		// query aligned to the reverse strand: the stored sequence is the reverse
		// complement of the original sequence
		// thus, reverse complement the stored sequence
		for (int i = 0; i < n; ++i) {
			cout << nuc_to_char(nuc_complement(bam_seqi(bam_get_seq(b), n - i - 1)));
		}
	} else {
		for (int i = 0; i < n; ++i) {
			cout << nuc_to_char(bam_seqi(bam_get_seq(b), i));
		}
	}
}

/*
* Print the phred-scaled quality for each nucleotide in the given BAM record.
* 
* @param b Pointer to BAM record containing qualities to be printed.
* @param original Whether the original (non reversed)
*	quality string should be returned for a reverse strand read. [false]
* @return none (prints the qualities in the record to stdout)
*/
// Wrong character offset? Currently prints from 'B' = phred 0 to 'L' = phred 10
void print_qual(const bam1_t *b, bool original = false) {
	if (b == NULL) return;
	const uint8_t* qual = bam_get_qual(b);
	if (original && bam_is_rev(b)) {
		// query aligned to the reverse strand: the stored quality is reversed
		// thus, reverse the stored quality
		int n = b->core.l_qseq;
		for (int i = 0; i < n; ++i) {
			cout << char(qual[n - i - 1] + 33);
		}
	} else {
		for (int i = 0; i < b->core.l_qseq; ++i) {
			cout << char((*qual) + 33);
			++qual;
		}
	}
}

/*
* Print the query name and sequence 
* in the given BAM record in FASTA format.
*/
void print_query_fasta(const bam1_t *b) {
	if (b == NULL) return;
	cout << '>' << bam_get_qname(b) << endl;
	print_seq(b, true); cout << endl;
}

/*
* Print the query name and sequence 
* in the given BAM record in FASTQ format.
*/
void print_query_fastq(const bam1_t *b) {
	if (b == NULL) return;
	cout << '@' << bam_get_qname(b) << endl;
	print_seq(b, true); cout << endl;
	cout << '+' << endl;
	print_qual(b, true); cout << endl;
}

/**
* Print the query nucleotide at the given reference position in the given BAM header,
* along with associated query name and position information.
* 
* @param b Pointer to BAM record containing the query nucleotide of interest
* @param pos Reference position of the query nucleotide
* @return none (Prints information about the selected query nucleotide to stdout)
*/
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

/**
* Print the names of the reference sequences listed in the
* BAM header stored in the given hts::fetcher object, up to a given index.
*
* @param f Properly initialized hts::fetcher containing the BAM header of interest.
* @param n Number of reference sequences to list, starting from the beginning.
* @return none (Prints names of reference sequences to stdout)
*/
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

/**
* Print a given number of reads from a BAM file stored in an hts::fetcher,
* along with their query name (read name), phred quality,
* cigar ops, and any auxiliary tags and their values.
* 
* @param f A properly initialized hts::fetcher containing the BAM data of interest.
* @param n The number of reads t to print
* @return none (Prints reads & read information to stdout)
*/

void print_n_alignment(fetcher& f, size_t n) {
	bam1_t *b = bam_init1();

	cout << "qname\tseq\tqual\tcigar\taux" << endl;

	for (size_t r = 0; r < n; ++r) {

		if (bam_read1(f.hf->fp.bgzf, b) > 0) {
			cout << bam_get_qname(b) << '\t';

			for (int i = 0; i < b->core.l_qseq; ++i) {
				cout << nuc_to_char(bam_seqi(bam_get_seq(b), i));
			}
			cout << '\t';

			const uint8_t* qual = bam_get_qual(b);
			for (int i = 0; i < b->core.l_qseq; ++i) {
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

