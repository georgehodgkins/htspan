#ifndef _HTSPAN_SNV_HPP_
#define _HTSPAN_SNV_HPP_

#include <fstream>
#include <sstream>
#include <stdexcept>

#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/khash.h"

#include "htspan/nucleotide.hpp"
#include "htspan/bam.hpp"

namespace hts {

using namespace std;

namespace snv {

// Initializes 'vdict' hash table alias used by HTSlib to store contig info
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)

struct record {
	// Human-readable name of reference sequence
	string chrom;
	// Read reference position
	int32_t pos;
	// Reference nucleotide
	nuc_t nt_ref;
	// Alternate nucleotide
	nuc_t nt_alt;
	/**
	* Non-fatal error codes:
	* 0=OK
	* 1=could not find RID in BAM header (or VCF header, for VCFs)
	* 2=could not unpack VCF 
	* 3=multiple-nucleotide variant in VCF
	* -1=Externally set error (record was read correctly)
	*/
	// Note that since vcf_read() returns the same value for EOF and other read errors, we cannot easily catch errors there
	int32_t err;
};

struct reader {

	virtual bool next (record &r) = 0;

};

struct tsv_reader : reader {
	ifstream f;

	string line;

	tsv_reader(const char* path) : f(path) {
	 	if (!f.is_open()) {
	 		throw runtime_error("Error: could not open input TSV file.");
	 	}

	 	// discard header line
	 	getline(f, line);
	 }

	/**
	 * Get next record.
	 */
	bool next(record& r) {
		// get next valid line
		do {
			if (f.eof()) return false;
			getline(f, line);
		} while (line.empty() || line[0] == '#');

		// process line
		istringstream ss(line);
		char char_ref, char_alt;
		ss >> r.chrom >> r.pos >> char_ref >> char_alt;
		r.err = 0;
		// convert from 1-based to 0-based
		r.pos -= 1;
		r.nt_ref = char_to_nuc(char_ref);
		r.nt_alt = char_to_nuc(char_alt);

		return true;
	}

	~tsv_reader() {
		f.close();
	}
};

struct vcf_reader : reader {
	// pointer to main VCF/BCF file object
	htsFile *hf;
	// pointer to VCF/BCF header object
	bcf_hdr_t *hdr;
	// pointer to VCF record buffer to read to/from
	bcf1_t *v;

	// the second parameter is for compatibility with the above's constructor signature, not used
	vcf_reader(const char* path) {
		// open HTS file handle
		hf = hts_open(path, "r");
		if (!hf) {
			throw runtime_error("Could not open input VCF file.");
		}

		// read VCF/BCF header from file
		hdr = bcf_hdr_read(hf);
		if (!hdr) {
			throw runtime_error("Could not read header from VCF file.");
		}

		// allocate space for VCF record
		v = bcf_init();
	}

	bool next(record &r) {
		// read in the record
		int status = bcf_read1(hf, hdr, v);
		if (status == -1) {// EOF or fatal reading error
			return false;
		}
		// unpack the record up to (including) the ALT field
		status = bcf_unpack(v, BCF_UN_STR);
		if (status < 0) {
			r.err = 2;
			return true;
		}
		// This class only handles /S/ NVs
		if (v->rlen > 1) {
			r.err = 3;
			return true;
		}
		// fill record fields
		r.err = 0;
		int rid = v->rid;
		// HTSlib automatically converts read chrom IDs from a VCF into int IDs,
		// even if that contig is not found in the header (it creates a dummy contig entry if that is the case)
		// We need the original chrom name to match with the BAM, so we have to go back in and find it
		//
		// find contig (chrom) string in the hash table
		typedef khash_t(vdict) vdict_t;
		vdict_t *contig_dict = (vdict_t*) hdr->dict[BCF_DT_CTG];
		r.chrom.clear();
		for (khint_t iter = kh_begin(contig_dict); iter != kh_end(contig_dict); ++iter) {
			bcf_idinfo_t contig = kh_val(contig_dict, iter);// val is a HTSlib struct
			if (contig.id == rid) {
				r.chrom = kh_key(contig_dict, iter);// key is c-string chrom name
				break;
			}
		}
		// if contig was not found
		if(r.chrom.empty()) {
			r.chrom = "VCF-unset";
		}
		r.pos = v->pos;// HTSlib internally converts from 1-based to 0-based
		r.nt_ref = char_to_nuc(v->d.allele[0][0]);
		r.nt_alt = char_to_nuc(v->d.allele[1][0]);
		return true;
	}

	~vcf_reader() {
		hts_close(hf);
		bcf_destroy1(v);
		bcf_hdr_destroy(hdr);
	}
};

}  // namespace snv

}  // namepsace hts

#endif  // _HTSPAN_SNV_HPP_
