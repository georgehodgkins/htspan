#ifndef _HTSPAN_SNV_READER_HPP_
#define _HTSPAN_SNV_READER_HPP_

#include <fstream>
#include <sstream>
#include <stdexcept>

#include "htslib/hts.h"
#include "htslib/vcf.h"

#include "htspan/nucleotide.hpp"
#include "htspan/bam.hpp"

//TODO: Documentation

namespace hts {

using namespace std;

namespace snv {

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

/**
* This reader reads from a TSV file containing the following columns:
* chrom	pos	ref	alt
* It must contain a header line.
*/
struct tsv_reader : reader {
	ifstream f;

	string line;

	record cached;

	tsv_reader(const char* path) {
	 	open(path);
	}

	/**
	* Opens the given path for reading as a TSV file.
	* Throws exception if the path cannot be opened.
	*/
	void open (const char* path) {
	 	f.open(path);

	 	if (!f.is_open()) {
	 		throw runtime_error("Error: could not open input TSV file.");
	 	}

	 	// discard header line
	 	// TODO: check header line
	 	getline(f, line);
	}

	/**
	 * Read next record and store it in the passed reference.
	 * A copy will also be cached internally until the next call to next(),
	 * accessible using get_underlying().
	 */
	bool next(record& r) {
		// get next valid line
		do {
			if (f.eof()) return false;
			getline(f, line);
		} while (line.empty() || line[0] == '#');
		// TODO: check for malformed lines

		// process line
		istringstream ss(line);
		char char_ref, char_alt;
		ss >> r.chrom >> r.pos >> char_ref >> char_alt;
		r.err = 0;
		// convert from 1-based to 0-based
		r.pos -= 1;
		r.nt_ref = char_to_nuc(char_ref);
		r.nt_alt = char_to_nuc(char_alt);

		cached = r;

		return true;
	}

	/**
	* Return cached record (most recently read).
	*/
	record get_underlying () const {
		return cached;
	}

	// Closes the attached file handle
	void close() {
		f.close();
	}

	// Destructor alias for close()
	~tsv_reader() {
		close();
	}
};

/**
* This reader uses HTSlib to read SNV data from VCF or BCF format files,
* optionally gzip or bgzip compressed. It will ignore all non-single-nucleotide variants.
*/
struct vcf_reader : reader {
	// pointer to main VCF/BCF file object
	htsFile *hf;
	// pointer to VCF/BCF header object
	bcf_hdr_t *hdr;
	// pointer to VCF record buffer to read to/from
	bcf1_t *v;

	// constructor alias for open()
	vcf_reader(const char* path) {
		open(path);
	}

	/**
	* Opens the given path for reading as a VCF file.
	* Throws exception if the file cannot be opened
	* or does not appear to be in VCF format.
	*/
	void open (const char* path) {
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

	/**
	* Read next record and store data from it in the passed reference.
	* A cached copy of the original bcf1_t object, which contains more fields
	* than a snv::record, will be available until the next call to next(),
	* accessible using get_underlying().
	*/
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
		r.chrom = bcf_hdr_id2name(hdr, v->rid);
		r.pos = v->pos;// HTSlib internally converts from 1-based to 0-based
		r.nt_ref = char_to_nuc(v->d.allele[0][0]);
		r.nt_alt = char_to_nuc(v->d.allele[1][0]);
		return true;
	}

	// Returns the last read bcf1_t object
	bcf1_t* get_underlying () const {
		return v;
	}

	/**
	* Closes the attached file handle and 
	* deallocates all memory allocated in open(),
	* including the object cache.
	*/
	void close() {
		if (hf != NULL) {
			hts_close(hf);
			hf = NULL;
		}
		if (v != NULL) {
			bcf_destroy1(v);
			v = NULL;
		}
		if (hdr != NULL) {
			bcf_hdr_destroy(hdr);
			hdr = NULL;
		}
	}

	// destructor alias for close()
	~vcf_reader() {
		close();
	}
};

}  // namespace snv

}  // namepsace hts

#endif  // _HTSPAN_SNV_READER_HPP_
