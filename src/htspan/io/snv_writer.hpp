#ifndef _HTSPAN_SNV_WRITER_HPP_
#define _HTSPAN_SNV_WRITER_HPP_

#include <stdexcept>
#include <cstring>
#include <fstream>
#include <queue>

#include "htslib/vcf.h"
#include "htslib/hts.h"

#include "snv_reader.hpp"
#include "../base_orient_bias_filter.hpp"

namespace hts {

namespace snv {

using namespace std;

// TODO: add base class for dummy methods

/**
* This snv_writer is used to write SNV records to a TSV file.
*/
struct tsv_writer {
	// underlying file stream
	ofstream f;

	// Constructor alias for open()
	// the to_copy parameter is a dummy for now,
	// to keep compatibility with the vcf_writer constructor
	// It is possible it will eventually be used if more features are added
	tsv_writer(const char* path, tsv_reader &to_copy) {
		open(path);
	}

	// Open the file at the given path for writing
	// NB: opening the file will destroy its former contents, if any
	void open (const char* path) {
		f.open(path, ios::trunc);
		if (!f.is_open()) {
			throw runtime_error("Error opening TSV file for output.");
		}
		// write header line
		f << "chrom\tpos\tref\talt" << endl;
	}

	// Write one record to the file (direct, immediate write)
	void write (record r, const char* filter_tag = NULL) {
		f << r.chrom << '\t' << r.pos+1 << '\t' << nuc_to_char(r.nt_ref) << '\t' << nuc_to_char(r.nt_alt) << endl;
	}

	// Close the underlying file stream
	void close() {
		f.close();
	}

	// This function does nothing, here for compatibilty
	// Will be moved to the base class when it is created
	void add_filter_tag (base_orient_bias_filter_f& filter) {}

	// Destructor alias for close()
	~tsv_writer() {
		close();
	}
};

/*
* This snv_writer is used to write SNV records to a VCF file (or related).
* Note that unlike the TSV writer, this class caches the header and written records
* and only writes them on a call to flush().
*/
struct vcf_writer {
	// Pointer to VCF file
	htsFile *hf;
	// Pointer to BCF header struct for output file
	bcf_hdr_t *hdr;
	// Linked header (to synchronize with
	bcf_hdr_t *linked;
	// FIFO of records to write
	queue<bcf1_t*> write_queue;//TODO: add pre-allocation from template header size

	// constructor alias for open()
	vcf_writer(const char* path, vcf_reader &to_copy) {
		open(path, to_copy.hf, to_copy.hdr);
	}

	// Opens a the given path for writing with the same format and compression mode as that
	// in the given htsFile; also links to the given vcf_header (needed to define contigs) 
	void open(const char* path, htsFile *f_templ, bcf_hdr_t *h_templ) {
		linked = h_templ;
		hdr = bcf_hdr_dup(linked);
		char mode[2];
		// write mode
		mode [0] = 'w';
		// determine compression from the file pointer directly,
		// since the API does not provide a method to do so
		switch(f_templ->format.compression) {
		case no_compression:
		default:
			mode[1] = 'u';
			break;
		case gzip:
			mode[1] = 'g';
			break;
		case bgzf:
			mode[1] = 'z';
		}
		hf = hts_open_format(path, mode, &f_templ->format);
		if (!hf) {
			throw runtime_error("Error opening VCF file for output.");
		}
	}

	// cache one VCF record for writing, adding a filter tag to it if requested
	void write(bcf1_t* rec, const char* filter_tag = NULL) {
		sync_header();// have to sync in case new contigs were added by a read
		if (filter_tag != NULL) {
			int filter_id = bcf_hdr_id2int(hdr, BCF_DT_ID, filter_tag);
			if (filter_id < 0) {
				throw runtime_error("Filter tag not found in header! Add the filter to the header first.");
			}
			int status = bcf_add_filter(hdr, rec, filter_id);
			if (status < 0) {
				throw runtime_error("Could not add filter tag to selected record.");
			}
		}
		// allocate a copy of the record (since the reader does not preserve records) and store a pointer to it
		write_queue.push(bcf_dup(rec));
	}

	// updates the writeout header with new contigs added to the linked header
	void sync_header() {
		bcf_hdr_merge(hdr, linked);
	}

	// Writes the header and all cached records to the file.
	// This should only be called when closing the file
	void flush() {
		sync_header();
		int status = bcf_hdr_write(hf, hdr);
		if (status != 0) {
			throw runtime_error("Error writing header to VCF file.");
		}
		while (!write_queue.empty()) {
			bcf1_t* twr = write_queue.front();
			status = bcf_write(hf, hdr, twr);
			if (status < 0) {
				throw runtime_error("Error writing record to VCF file.");
			}
			write_queue.pop();
			// deallocate copy of record
			bcf_destroy(twr);
		}
	}

	// Adds the tag from the given filter to the file header
	void add_filter_tag (base_orient_bias_filter_f& filter) {
		sync_header();
		int status = bcf_hdr_printf(hdr, "##FILTER=<ID=%s,Description=\"%s\">", filter.text_id, filter.get_description());
		if (status < 0) {
			throw runtime_error("Could not add filter tag to VCF header.");
		}
	}

	// flushes the file and then closes it
	void close() {
		if (hf != NULL) {
			flush();
			linked = NULL;
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