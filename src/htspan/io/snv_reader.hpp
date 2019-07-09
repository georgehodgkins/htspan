#ifndef _HTSPAN_SNV_READER_HPP_
#define _HTSPAN_SNV_READER_HPP_

#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <ostream>

#include "htslib/hts.h"
#include "htslib/vcf.h"

#include "htspan/nucleotide.hpp"
#include "htspan/bam.hpp"

#include "snv.hpp"
#include "frontend/cstring.hpp"

// TODO: Make readers uniformly follow the path of 1) validate data 2) store in cache 3) copy to output record

namespace hts {

using namespace std;

namespace snv {

/*
* Base class for SNV readers; exposes the methods used in calling code.
*/
struct reader {

	/**
	* Read the next record from the file and 
	* fill the corresponding fields of the passed
	* snv::record.
	*
	* Returns true if there is more data, and false otherwise.
	* The error() method indicates whether an error occurred during
	* the most recent read (see below); in addition, the record will be
	* set to null (test with is_null()) if an error occurred. If EOF was reached,
	* the record will /not/ be set to null.
	*/
	virtual bool next (record &r) = 0;

	int err;

	size_t S;

	/**
	* Get the error code from the last read.
	* Codes:
	* -2: Reader is unininitialized (has not read a record since being opened)
	* -1: Externally set error (no problem with read)
	*  0: No error
	*  1: Read alternate allele was invalid or zero-length (deletion)
	*  2: Could not unpack BCF record
	*  3: Read alternate allele was greater than one nucleotide in length
	*/
	int error () const {
		return err;
	}

	reader () {
		err = -2;
		S = 0;
	}

	/**
	* Returns the format of the read file as a 
	* FMTFLAGS_T enum (defined in snv.hpp).
	*/
	virtual FMTFLAGS_T get_format () const = 0;

	/**
	* Opens the file at the given path for reading.
	* The VCF reader will throw an exception if attempting to open an
	* incompatible format, but the TSV reader will not.
	*/
	virtual void open (const char* path) = 0;

	/**
	* Returns a const reference to the
	* cached copy of the last read record.
	*/
	virtual const record& get_cached () const = 0;

	/**
	* Closes the attached file handle.
	*/
	virtual void close () = 0;

	virtual ~reader() {}

};

/**
* This reader reads from a TSV file containing the following columns:
* chrom	pos	ref	alt
* It must contain a header line.
*/
struct tsv_reader : reader {
	// underlying file stream
	ifstream f;
	// last line read from stream
	string line;
	// cached last-read record
	record cached;
	// index to look for the next alt at, for multiple alts
	size_t alt_i;
	// max index an alt can be present at
	size_t alt_max;

	tsv_reader(const char* path) : reader() {
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
	 *
	 * Return value indicates whether the EOF was reached.
	 * Errors return true and set the error flag instead.
	 */
	bool next(record& rec) {
		// if not at EOF, more alt nucs to read
		if (!err && alt_i < alt_max) {
			char nuc_alt = line[alt_i];
			// check for multi-nuc alt
			if (alt_i < alt_max - 1 && line[alt_i + 1] != ',' && line[alt_i + 1] != '\t') {
				err = 3;
				return true;
			}
			cached.nt_alt = char_to_nuc(nuc_alt);
			// check for zero-nuc alt/invalid character
			if (!nuc_is_canonical(cached.nt_alt)) {
				err = 1;
				return false;
			}
			// skip the comma
			alt_i += 2;
			// copy data to return record
			rec = cached;
			return true;
		} else {// get next valid line
			err = 0;
			// skip empty lines and comments
			do {
				if (f.eof()) return false;
				getline(f, line);
			} while (line.empty() || line[0] == '#');
			// TODO: check for malformed lines
			// clear previous
			rec.clear();
			cached.clear();

			// read first three fields (alt nucs are handled separately)
			istringstream ss(line);
			char char_ref, char_alt;
			ss >> cached.chrom >> cached.pos >> char_ref;
			// ref length < 1
			if (char_ref == '-' || char_ref == '.') {
				err = 1;
				return true;
			}
			// ref length > 1
			if (ss.peek() != '\t') {
				err = 3;
				return true;
			}

			// convert from 1-based to 0-based
			cached.pos -= 1;
			cached.nt_ref = char_to_nuc(char_ref);

			// account for multiple alts:
			// find the last tab
			// check if the info field is present

			alt_i = line.find_last_of("\t") + 1;
			size_t n_tabs = str_count_char(line.c_str(), '\t');
			if (n_tabs > 3) { // an info field is present, so find the next tab back
				alt_max = alt_i;
				alt_i = line.find_last_of("\t", alt_max - 2) + 1;
			} else { // no info field, so max is EOL
				alt_max = line.size();
			}
			char_alt = line[alt_i];
			cached.nt_alt = char_to_nuc(char_alt);

			// alt length < 1 (or invalid character)
			if (!nuc_is_canonical(cached.nt_alt)) {
				err = 1;
				return true;
			}
			// alt length > 1
			if (alt_i < alt_max - 1 && line[alt_i + 1] != ',' && line[alt_i + 1] != '\t') {
				err = 3;
				return true;
			}
			// if multiple nuc variants are excluded, we can
			// assume the next alt nuc would be the first + 2 (skip the comma) 
			alt_i += 2;

			// copy read data to returned object
			rec = cached;
		}
		return true;
	}

	FMTFLAGS_T get_format () const {
		return F_TSV;
	}

	/**
	* Return a reference to cached record (most recently read).
	*/
	const record& get_cached () const {
		return cached;
	}

	// Closes the attached file handle
	void close() {
		f.close();
		err = -2;
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
	// cached record (last read)
	// the internal bcf1_t record serves as the buffer for reading
	record cached;
	// counter for reading multiple alleles from a single record
	uint16_t rd_als;

	// constructor alias for open()
	vcf_reader(const char* path) : reader() {
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
		cached.v = bcf_init();
	}

	/**
	* Read next record and store data from it in the passed reference.
	* A cached copy of the original bcf1_t object, which contains more fields
	* than a snv::record, will be available until the next call to next(),
	* accessible using get_underlying().
	*/
	bool next(record &r) {
		// still alt alleles to read from the current record
		if (!err && rd_als < cached.v->n_allele) {
			// clear the record only if it contains different 
			if (r.v == NULL || r.v->pos != cached.v->pos || r.v->rid != cached.v->rid) {
				r.clear();
			}
			size_t alen = strlen(cached.v->d.allele[rd_als]);
			
			// zero-nuc alt
			if (alen < 1) {
				err = 1;
				r.clear();
				return true;
			}

			// multi-nuc alt
			if (alen > 1) {
				err = 3;
				r.clear();
				return true;
			}

			// if checks pass, fill record fields
			cached.chrom = bcf_hdr_id2name(hdr, cached.v->rid);
			cached.pos = cached.v->pos;// HTSlib internally converts from 1-based to 0-based
			cached.nt_ref = char_to_nuc(cached.v->d.allele[0][0]);
			cached.nt_alt = char_to_nuc(cached.v->d.allele[rd_als][0]);
			r = cached;
			++rd_als;
			return true;

		} else { // read a new record
			rd_als = 0;
			err = 0;

			// read in the record
			int status = bcf_read1(hf, hdr, cached.v);
			if (status == -1) {// EOF or other fatal reading error
				return false;
			}

			r.clear();

			// unpack the record up to (including) the ALT field
			status = bcf_unpack(cached.v, BCF_UN_STR);
			if (status < 0) {
				err = 2;
				return true;
			}

			// Multi-nuc ref
			if (cached.v->rlen > 1) {
				err = 3;
				return true;
			}

			// fill record fields
			cached.chrom = bcf_hdr_id2name(hdr, cached.v->rid);
			cached.pos = cached.v->pos;// HTSlib internally converts from 1-based to 0-based
			cached.nt_ref = char_to_nuc(cached.v->d.allele[0][0]);

			// catch zero-length ref and misc invalid chars
			if (!nuc_is_canonical(cached.nt_ref)) {
				err = 1;
				return true;
			}

			// length of alt
			size_t alen = strlen(cached.v->d.allele[1]);

			// multi-nuc alt
			if (alen > 1) {
				err = 3;
				return true;
			}

			cached.nt_alt = char_to_nuc(cached.v->d.allele[1][0]);

			// check for zero-length nuc or invalid char
			if (!nuc_is_canonical(cached.nt_alt)) {
				err = 1;
				return false;
			}
			rd_als = 2;

			// copy data to returned object
			r = cached;
			
			return true;
		}
	}

	FMTFLAGS_T get_format () const {
		if (hf->is_bgzf) {
			return static_cast<FMTFLAGS_T>(F_BGZF | F_VCF);
		} else {
			return F_VCF;
		}
	}

	// Returns a reference to the last read bcf1_t object
	const record& get_cached () const {
		return cached;
	}

	/**
	* Closes the attached file handle and 
	* deallocates all memory allocated in open().
	*
	* Note that snv::record has a destructor and so
	* `cached` does not need to be destroyed explicitly.
	*/
	void close() {
		if (hf != NULL) {
			hts_close(hf);
			hf = NULL;
		}
		if (hdr != NULL) {
			bcf_hdr_destroy(hdr);
			hdr = NULL;
		}
		err = -2;
	}

	// destructor alias for close()
	~vcf_reader() {
		close();
	}
};

}  // namespace snv

}  // namepsace hts

#endif  // _HTSPAN_SNV_READER_HPP_
