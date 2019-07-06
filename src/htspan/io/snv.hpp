#ifndef _HTSPAN_SNV_HPP_
#define _HTSPAN_SNV_HPP_

#include "frontend/cstring.hpp"

namespace hts {

namespace snv {

enum FMTFLAGS_T {
	F_NULL = 0x00,
	F_TSV = 0x01,
	F_VCF = 0x02,
	F_BGZF = 0x04
};

FMTFLAGS_T operator| (FMTFLAGS_T a, FMTFLAGS_T b) {
	return static_cast<FMTFLAGS_T>(static_cast<int>(a) | static_cast<int>(b));
}

FMTFLAGS_T operator& (FMTFLAGS_T a, FMTFLAGS_T b) {
	return static_cast<FMTFLAGS_T>(static_cast<int>(a) & static_cast<int>(b));
}

/**
* Converts command-line option strings to a FMTFLAGS_T enum.
* Uses typearg if it is not null; otherwise, attempts to deduce
* type from filename extension, returns F_NULL on failure to deduce.
* 
* @param filename Filename passed on command line
* @param typearg Type argument passed on command line
* @return Appropriate format flag, or F_NULL for invalid inputs
*/
FMTFLAGS_T opts_to_fmt (const char* filename, const char* typearg) {
	if (typearg) {
		if (strcmpi(typearg, "tsv") == 0) {
			return F_TSV;
		} else if (strcmpi(typearg, "vcf") == 0) {
			return F_VCF;
		} else if (strcmpi(typearg, "vcf-bgz") == 0) {
			return (F_VCF | F_BGZF);
		} else { // invalid type argument
			return F_NULL;
		}
	} else if (filename) {
		const char* xtn = strrchr(filename, '.');
		if (xtn) {
			++xtn;
			if (strcmpi(xtn, "snv") == 0 ||
					strcmpi(xtn, "tsv") == 0) {
				return F_TSV;
			} else if (strcmpi(xtn, "vcf") == 0 ||
					strcmpi(xtn, "bcf") == 0) {
				return F_VCF;
			} else if (strcmpi(xtn, "gz") ||
					strcmpi(xtn, "bgz")) {
				return (F_VCF | F_BGZF);
			} else { // unrecognized extension on filename
				return F_NULL;
			}
		} else { // no '.' in filename
			return F_NULL;
		}
	} else { // null extension and filename
		return F_NULL;
	}
}

/**
* Return the file extension corresponding to the
* passed FMTFLAGS_T enum.
*/
const char* fmt_to_xtn (FMTFLAGS_T fmt) {
	if (fmt & F_VCF) {
		if (fmt & F_BGZF) {
			return ".vcf.gz";
		} else {
			return ".vcf";
		}
	} else if (fmt & F_TSV) {
		return ".snv";
	} else { // F_NULL
		return "";
	}
}

/**
* This class represents a single-nucleotide variant,
* storing its human-readable name (usually a chromosome number),
* reference position, and reference and alternative nucleotides.
* If the record was read from a VCF/BCF file, it will also contain a pointer
* to a copy of the HTSlib bcf1_t datatype, which contains much more information
* about the variant.
*/
struct record {
	// Human-readable name of reference sequence
	string chrom;
	// Read reference position
	int32_t pos;
	// Reference nucleotide
	nuc_t nt_ref;
	// Alternate nucleotide
	nuc_t nt_alt;
	// pointer to a corresponding BCF record (NULL if not present)
	// NB: must be unpacked with bcf_unpack for most data to be accesible
	bcf1_t *v;

	// parameter constructor ( v must be set directly)
	record (const char* c, int32_t p, char r, char a)
		: chrom(c), pos(p), nt_ref(char_to_nuc(r)), nt_alt(char_to_nuc(a)) {
			v = NULL;
		}

	// constructs a null instance of the class (calls clear())
	record () {
		v = NULL; // must be set here to indicate it is not allocated for clear()
		clear();
	}

	// copy constructor (calls operator=())
	// TODO: this should really be the other way around
	record (const record &r) {
		v = NULL;
		operator=(r);
	}	

	// Prints a human-readable string describing the SNV
	string to_string () const {
		if (is_null()) {
			return "EMPTY";
		}
		ostringstream ss;
		ss << nuc_to_char(nt_ref) << '>' << nuc_to_char(nt_alt) << " @ " << chrom << ':' << pos;
		if (v != NULL) {
			ss << " [vcf]";
		}
		return ss.str();
	}

	/*
	* Nulls the object, setting all fields to dummy values.
	* If bool true is passed as an argument, it will not
	* deallocate the internal bcf record (but will still
	* clear its data).
	*/
	void clear (bool preserve_bcf = false) {
		chrom = "";
		pos = -1337;
		nt_ref = nuc_N;
		nt_alt = nuc_N;
		if (v != NULL) {
			if (preserve_bcf) {
				bcf_clear(v);
			} else {
				bcf_destroy1(v);
				v = NULL;
			}
		}
	}

	/*
	* Deep copy operator for the class.
	*/
	void operator= (const record &p) {
		if (p.v == NULL) {
			clear();
		} else if (v == NULL) {
			clear();
			v = bcf_dup(p.v);
		} else {
			clear(true); // clears but does not deallocate bcf record
			bcf_copy(v, p.v);
		}
		chrom = p.chrom;
		pos = p.pos;
		nt_ref = p.nt_ref;
		nt_alt = p.nt_alt;
	}

	// Returns whether the object is null.
	bool is_null () const {
		// default constructor constructs the null case
		return operator==( record() );
	}

	~record () {
		clear();
	}

	// Comparison overloads. Note that they do not consider the internal BCF record.
	bool operator== (const record &p) const {
		return chrom == p.chrom && pos == p.pos && nt_ref == p.nt_ref && nt_alt == p.nt_alt;
	}

	bool operator!= (const record &p) const {
		return !(operator==(p));
	}

};

} // namespace snv

} // namespace hts

#endif // _HTSPAN_SNV_HPP_