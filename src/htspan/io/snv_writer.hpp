#ifndef _HTSPAN_SNV_WRITER_HPP_
#define _HTSPAN_SNV_WRITER_HPP_

#include <stdexcept>
#include <cstring>
#include <fstream>
#include <queue>

#include "htslib/vcf.h"
#include "htslib/hts.h"

#include "../base_orient_bias_filter.hpp"
#include "snv.hpp"

namespace hts {

namespace snv {

using namespace std;

struct writer {

	virtual void add_filter_tag(const base_orient_bias_filter_f &filter) {};

	virtual void write(const record &rec) {
		if (rec.is_null()) {
			throw invalid_argument("Invalid SNV record passed to writer! Aborting.");
		}
	}

	virtual void write_filter_failed (const record &r, const base_orient_bias_filter_f &filter) = 0;

	virtual void open(const char* path, const FMTFLAGS_T fmt) = 0;

	virtual void close() = 0;

};

/**
* This snv_writer is used to write SNV records to a TSV file.
*/
struct tsv_writer : writer {
	// underlying file stream
	ofstream f;

	// Constructor alias for open()
	// the to_copy parameter is a dummy for now,
	// to keep compatibility with the vcf_writer constructor
	// It is possible it will eventually be used if more features are added
	tsv_writer(const char* path, const FMTFLAGS_T fmt) {
		open(path, fmt);
	}

	// Open the file at the given path for writing
	// NB: opening the file will destroy its former contents, if any
	void open (const char* path, const FMTFLAGS_T fmt) {
		if (!(fmt & F_SNV)) {// if F_SNV is not set
			throw runtime_error("TSV writer format must be F_SNV.");
		}
		f.open(path, ios::trunc);
		if (!f.is_open()) {
			throw runtime_error("Error opening TSV file for output.");
		}
		// write header line
		f << "chrom\tpos\tref\talt" << endl;
	}

	// Write one record to the file (direct, immediate write)
	void write (const record &r) {
		writer::write(r);
		f << r.chrom << '\t' << r.pos+1 << '\t' << nuc_to_char(r.nt_ref) << '\t' << nuc_to_char(r.nt_alt) << endl;
	}

	// in this class, this method intentionally does nothing,
	// since record should not be written back if filter fails
	void write_filter_failed (const record &r, const base_orient_bias_filter_f &filter) {}

	// Close the underlying file stream
	void close() {
		f.close();
	}

	// Destructor alias for close()
	~tsv_writer() {
		close();
	}
};

/*
* This snv_writer is used to write SNV records to a VCF file (or related).
* Note that unlike the TSV writer, this class caches the header and written records
* and only writes them on a call to flush().
* TODO: change to immediate write
*
* The class must be linked to a header read from an existing BCF file, 
* from which it copies file format data and header fields.
* TODO: make this optional
*/
struct vcf_writer : writer {
	// Pointer to VCF file
	htsFile *hf;
	// Pointer to BCF header struct for output file
	bcf_hdr_t *hdr;
	// FIFO of records to write
	queue<bcf1_t*> write_queue;//TODO: add pre-allocation from template header size

	// constructor alias for open()
	vcf_writer(const char* path, FMTFLAGS_T fmt) {
		open(path, fmt);
	}

	// Opens a the given path for writing with the same format and compression mode as that
	// in the given htsFile; also links to the given vcf_header (needed to define contigs) 
	void open(const char* path, FMTFLAGS_T fmt) {
		if (!(fmt & F_VCF)) {
			throw runtime_error("VCF writer format must be F_VCF.");
		}
		char mode[3];
		// write mode
		mode [0] = 'w';
		// determine compression from the file pointer directly,
		// since the API does not provide a method to do so
		if (fmt & F_BGZF) {
			mode[1] = 'z';
		}
		else {
			mode[1] = 'u';
		}
		mode[3] = '\0';
		hdr = bcf_hdr_init("w");
		if (!hdr) {
			throw runtime_error("Error initializing VCF header.");
		}
		hf = hts_open(path, mode);
		if (!hf) {
			throw runtime_error("Error opening VCF file for output.");
		}
	}

	// cache one VCF record for writing after adding a filter tag to it
	void write_filter_failed(const record &rec, const base_orient_bias_filter_f &filter) {
		int filter_id = bcf_hdr_id2int(hdr, BCF_DT_ID, filter.text_id);
		if (filter_id < 0) {
			throw runtime_error("Filter tag not found in header! Add the filter to the header first with add_filter_tag.");
		}
		int status = bcf_add_filter(hdr, rec.v, filter_id);
		if (status < 0) {
			throw runtime_error("Could not add filter tag to selected record.");
		}
		// cache record for writing
		write(rec);
	}

	void write(const record &rec) {
		writer::write(rec);
		bcf1_t *twr = bcf_dup(rec.v);
		// check that a corresponding contig exists
		if (!hdr->id[BCF_DT_CTG] || bcf_hdr_name2id(hdr, rec.chrom.c_str()) != twr->rid) {
			// if it does not exist, add it to the dictionary
			bcf_hdr_printf(hdr, "##contig=<ID=%s>", rec.chrom.c_str());
			twr->rid = bcf_hdr_name2id(hdr, rec.chrom.c_str());
		}
		write_queue.push(bcf_dup(twr));
	}

	// Writes the header and all cached records to the file.
	// This should only be called when closing the file
	void flush() {
		int status = bcf_hdr_write(hf, hdr);
		if (status != 0) {
			throw runtime_error("Error writing header to VCF file.");
		}
		while (!write_queue.empty()) {
			bcf1_t *twr = write_queue.front();
			status = bcf_write(hf, hdr, twr);
			if (status < 0) {
				throw runtime_error("Error writing record to VCF file.");
			}
			write_queue.pop();
			// deallocate copy of record just written
			bcf_destroy(twr);
		}
	}

	void copy_header (const bcf_hdr_t *src) {
		// adds header lines from source header
		bcf_hdr_merge(hdr, src);
		// copy over samples manually
		for (size_t n = 0; n < bcf_hdr_nsamples(src); ++n) {
			bcf_hdr_add_sample(hdr, src->samples[n]);
		}
	}

	// Adds the tag from the given filter to the file header
	void add_filter_tag (base_orient_bias_filter_f &filter) {
		int status = bcf_hdr_printf(hdr, "##FILTER=<ID=%s,Description=\"%s\">", filter.text_id, filter.get_description());
		if (status < 0) {
			throw runtime_error("Could not add filter tag to VCF header.");
		}
	}

	// flushes the file and then closes it
	void close() {
		if (hf != NULL) {
			flush();
			hts_close(hf);
			hf = NULL;
			bcf_hdr_destroy(hdr);
			hdr = NULL;
		}
	}

	// destructor alias for close()
	~vcf_writer () {
		close();
	}
};

}// namespace snv

}// namespace hts

#endif // _HTSPAN_SNV_WRITER_HPP_