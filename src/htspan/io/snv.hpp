#ifndef _HTSPAN_SNV_HPP_
#define _HTSPAN_SNV_HPP_

namespace hts {

namespace snv {

enum FMTFLAGS_T {
	F_NULL = 0x00,
	F_SNV = 0x01,//TODO: change to F_TSV
	F_VCF = 0x02,
	F_BGZF = 0x04
};

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

	// Prints a human-readable string describing the SNV
	string to_string () const {
		if (is_null()) {
			return "EMPTY";
		}
		ostringstream ss;
		ss << "Chrom: " << chrom << " Pos: " << pos << " Ref: " << nuc_to_char(nt_ref) << " Alt: " << nuc_to_char(nt_alt);
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
	// TODO: fix allocation in this method
	void operator= (const record &p) {
		if (p.v == NULL) {
			clear();
		} else if (v == NULL) {
			clear();
			v = bcf_dup(p.v);
		} else {
			// do not deallocate bcf record unless it differs
			clear(true);
			if (p.v->rid != v->rid || p.v->pos != v->pos) {
				bcf_destroy(v);
				v = bcf_dup(p.v);
			}
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