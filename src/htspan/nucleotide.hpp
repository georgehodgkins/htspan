#ifndef _HTSPAN_NUCLEOTIDE_HPP_
#define _HTSPAN_NUCLEOTIDE_HPP_

typedef uint8_t nuc_t;

/**
 * Extended nucleotide enum.
 */
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

/**
 * Whether nucleotide is A, C, G, or T.
 */
bool nuc_is_canonical(nuc_t x) {
	return x <= nuc_T;
}

/**
 * Complement nucleotide.
 */
nuc_t nuc_complement(nuc_t x) {
	switch (x) {
		case nuc_A: return nuc_T;
		case nuc_C: return nuc_G;
		case nuc_G: return nuc_C;
		case nuc_T: return nuc_A;
		default: return x;
	}
}

/**
 * Convert from nucleotide enum to ASCII char.
 */
char nuc_to_char(nuc_t x) {
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

/**
 * Convert from ASCII char to nucleotide enum.
 */
nuc_t char_to_nuc(char x) {
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

/**
 * Check whether two nucleotides are equal.
 * 
 * If either nucleotide is N, then they are equal.
 */
bool nuc_equal(nuc_t x, nuc_t y) {
	if (x == nuc_N || y == nuc_N) return true;
	if (x == y) return true;
	return false;
}

#endif  // _HTSPAN_NUCLEOTIDE_HPP_

