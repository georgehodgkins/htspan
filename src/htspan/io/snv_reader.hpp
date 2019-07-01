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
	// pointer to a corresponding BCF record (NULL if not present)
	// NB: must be unpacked with bcf_unpack for most data to be accesible
	bcf1_t *v;

	record (const char* c, int32_t p, char r, char a)
		: chrom(c), pos(p), nt_ref(char_to_nuc(r)), nt_alt(char_to_nuc(a)) {
			v = NULL;
		}

	record () {
		v = NULL; // must be set here to indicate it is not allocated for clear()
		clear();
	}	

	// used for debugging
	string to_string () const {
		ostringstream ss;
		ss << "Chrom: " << chrom << " Pos: " << pos << " Ref: " << nuc_to_char(nt_ref) << " Alt: " << nuc_to_char(nt_alt);
		return ss.str();
	}

	void clear (bool preserve_bcf = false) {
		chrom = "";
		pos = -1337;
		nt_ref = nuc_N;
		nt_alt = nuc_N;
		if (v != NULL && !preserve_bcf) {
			bcf_destroy1(v);
			v = NULL;
		}
	}

	void operator= (const record &p) {
		if (p.v == NULL) {
			clear();
		} else {
			clear(true);// do not deallocate bcf record
		}
		chrom = p.chrom;
		pos = p.pos;
		nt_ref = p.nt_ref;
		nt_alt = p.nt_alt;
		if (p.v != NULL) {
			v = bcf_dup(p.v);
		}
	}

	~record () {
		clear();
	}

	bool operator== (const record &p) const {
		return chrom == p.chrom && pos == p.pos && nt_ref == p.nt_ref && nt_alt == p.nt_alt;
	}

	bool operator!= (const record &p) const {
		return !(operator==(p));
	}

};


struct reader {

	virtual bool next (record &r) = 0;

	int err;

	int error () const {
		return err;
	}

};

/**
* This reader reads from a TSV file containing the following columns:
* chrom	pos	ref	alt
* It must contain a header line.
*/
struct tsv_reader : reader {
	ifstream f;

	string line;
	// cached last-read record
	record cached;
	// index to look for the next alt at, for multiple alts
	size_t alt_i;

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
	 	err = -2;
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
		// ss is a class member so its contents are preserved
		// if not at EOF, more alt nucs to read
		if (!err && alt_i < line.size()) {
			char nuc_alt = line[alt_i];
			// check for multi-nuc alt
			if (alt_i < line.size() - 1 && line[alt_i + 1] != ',') {
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
			rec.clear();
			cached.clear();
			// skip empty lines and comments
			do {
				if (f.eof()) return false;
				getline(f, line);
			} while (line.empty() || line[0] == '#');
			// TODO: check for malformed lines

			// process line
			// ss is a class member
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
			// index starts at the last whitespace char + 1,
			// assumed to be the first alt nuc
			alt_i = line.find_last_of("\t ") + 1;
			char_alt = line[alt_i];
			cached.nt_alt = char_to_nuc(char_alt);
			// alt length < 1 (or invalid character)
			if (!nuc_is_canonical(cached.nt_alt)) {
				err = 1;
				return true;
			}
			// alt length > 1
			if (alt_i < line.size() - 1 && line[alt_i + 1] != ',') {
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

	/**
	* Return a reference to cached record (most recently read).
	*/
	const record& get_cached () const {
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
	// cached record (last read)
	// the internal bcf1_t record serves as the buffer for reading
	record cached;
	// counter for reading multiple alleles from a single record
	uint16_t rd_als;

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
				return true;
			}

			// multi-nuc alt
			if (alen > 1) {
				err = 3;
				return true;
			}

			// if checks pass, fill record fields
			r.chrom = bcf_hdr_id2name(hdr, cached.v->rid);
			r.pos = cached.v->pos;// HTSlib internally converts from 1-based to 0-based
			r.nt_ref = char_to_nuc(cached.v->d.allele[0][0]);
			r.nt_alt = char_to_nuc(cached.v->d.allele[rd_als][0]);
			if (r.v == NULL) {
				r.v = bcf_dup(cached.v);
			}
			++rd_als;
			return true;

		} else { // read a new record
			r.clear();
			cached.clear(true); // do not deallocate buffer
			rd_als = 0;
			err = 0;

			// read in the record
			int status = bcf_read1(hf, hdr, cached.v);
			if (status == -1) {// EOF or other fatal reading error
				return false;
			}

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

	// Returns a reference to the last read bcf1_t object
	const record& get_cached () const {
		return cached;
	}

	/**
	* Closes the attached file handle and 
	* deallocates all memory allocated in open().
	*
	* Note that snv::record has a destructor and
	* does not need to be destroyed explicitly.
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
	}

	// destructor alias for close()
	~vcf_reader() {
		close();
	}
};

}  // namespace snv

}  // namepsace hts

#endif  // _HTSPAN_SNV_READER_HPP_
