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

// TODO: make multiple alt allele split write back to VCF correctly

/**
* Base class for SNV writers; exposes the methods used in calling code.
*/
struct writer {
	// Adds the tag from the passed filter or filter type to the tag list, if the reader supports annotation
	virtual void add_filter_tag(const base_orient_bias_filter_f &filter) {};

	//
	virtual void add_numeric_info(const char* name, const char* descr) = 0;

	// Writes a record to the file, or caches it for writing
	virtual void write(const record &rec, const char* info_name = NULL, const double info_val = -1.0) {
		if (rec.is_null()) {
			throw invalid_argument("Invalid SNV record passed to writer! Aborting.");
		}
	}

	// Writes a record, noting that it failed the filter
	// Depending on the behavior of the writer, this may not actually write the record
	virtual void write_filter_failed (const record &r, const base_orient_bias_filter_f &filter,
		const char* info_name = NULL, const double info_val = -1.0) = 0;

	// Open a file at the given path with the given format
	// FMTFLAGS_T is defined in snv.hpp
	virtual void open(const char* path, const FMTFLAGS_T fmt) = 0;

	// Close the attached file
	virtual void close() = 0;

	virtual ~writer() {}

};

/**
* This snv_writer is used to write SNV records to a TSV file.
*/
struct tsv_writer : writer {
	// underlying file stream
	ofstream f;
	// info field (only a single field is supported for this reader atm)
	string info_name;
	// header line, emptied once it is written
	string header_line;

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
		if (!(fmt & F_TSV)) {// if F_SNV is not set
			throw runtime_error("TSV writer format must be F_TSV.");
		}
		f.open(path, ios::trunc);
		if (!f.is_open()) {
			throw runtime_error("Error opening TSV file for output.");
		}
		// set header line
		header_line = "chrom\tpos\tref\talt\n";
	}

	void add_numeric_info (const char* name, const char* descr) {
		if (!info_name.empty()) {
			throw runtime_error("Cannot add more than one info field to this TSV reader!");
		}
		if (header_line.empty()) {
			throw runtime_error("Cannot add an info field after writing the header line!");
		}
		info_name = name;
		header_line.erase(header_line.end() - 1);//remove newline
		header_line += ('\t' + info_name + '\n');// append new field
	}

	void write_header_line () {
		f << header_line;
		header_line.clear();
	}

	// Write one record to the file (direct, immediate write)
	void write (const record &r, const char* info_field = NULL, const double info_val = -1.0) {
		writer::write(r);
		// write the header line if it has not yet been written
		if (!header_line.empty()) {
			write_header_line();
		}
		// write the line with an info field if it is given, without one otherwise
		if (info_field != NULL) {
			if (strcmp(info_field, info_name.c_str()) != 0) {
				throw runtime_error("Could not find given info field in the reader!");
			}
			f << r.chrom << '\t' << r.pos+1 << '\t' << nuc_to_char(r.nt_ref) << '\t' << nuc_to_char(r.nt_alt) << '\t' << info_val << endl;
		} else {
			f << r.chrom << '\t' << r.pos+1 << '\t' << nuc_to_char(r.nt_ref) << '\t' << nuc_to_char(r.nt_alt) << endl;
		}
	}

	// in this class, this method intentionally does nothing,
	// since record should not be written back if filter fails
	void write_filter_failed (const record &r, const base_orient_bias_filter_f &filter, 
		const char* info_name = NULL, const double info_val = -1.0) {}

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
*
* This is necessary because the header may change (added contigs and samples) when writing 
* records, and in order to write the header first in the file with HTSlib it must be the
* first thing written.
*
* The class must be linked to a header read from an existing BCF file, 
* from which it copies file format data and header fields.
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
		mode[2] = '\0';
		hdr = bcf_hdr_init("w");
		if (!hdr) {
			throw runtime_error("Error initializing VCF header.");
		}
		hf = hts_open(path, mode);
		if (!hf) {
			throw runtime_error("Error opening VCF file for output.");
		}
	}

	void add_numeric_info(const char* name, const char* descr) {
		bcf_hdr_printf(hdr, "##INFO=<ID=%s,Number=1,Type=Float,Description=\"%s\",Source=\"htspan-orient-bias\">", name, descr);
		bcf_hdr_sync(hdr); // add the newly added field to the ID list
	}

	// cache one VCF record for writing after adding a filter tag to it
	void write_filter_failed(const record &rec, const base_orient_bias_filter_f &filter,
			const char* info_name = NULL, const double info_val = -1.0) {
		int filter_id = bcf_hdr_id2int(hdr, BCF_DT_ID, filter.text_id);
		if (filter_id < 0) {
			throw runtime_error("Filter tag not found in header! Add the filter to the header first with add_filter_tag.");
		}
		int status = bcf_add_filter(hdr, rec.v, filter_id);
		if (status < 0) {
			throw runtime_error("Could not add filter tag to selected record.");
		}
		// cache record for writing
		write(rec, info_name, info_val);
	}

	// Caches a record for writing
	void write(const record &rec, const char* info_name = NULL, const double info_val = -1.0) {
		writer::write(rec);
		bcf1_t *twr = bcf_dup(rec.v); // this copy is destroyed after being written in flush()
		// check that a corresponding contig exists
		if (!hdr->id[BCF_DT_CTG] || bcf_hdr_name2id(hdr, rec.chrom.c_str()) != twr->rid) {
			// if it does not exist, add it to the dictionary
			bcf_hdr_printf(hdr, "##contig=<ID=%s>", rec.chrom.c_str());
			twr->rid = bcf_hdr_name2id(hdr, rec.chrom.c_str());
		}
		// add the info field if one is passed
		if (info_name != NULL) {
			if (bcf_hdr_id2int(hdr, BCF_DT_ID, info_name) == -1) {
				throw runtime_error("Info tag not found in header! Add the tag to the header first with add_numeric_info.");
			}
			bcf_update_info_float(hdr, twr, info_name, (void*) &info_val, 1);
		}
		// push the record onto the write queue
		write_queue.push(twr);
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

	/**
	* Copies the header and sample data from the given header.
	*
	* This is necessary if you are going to be copying
	* records from an existing file, to avoid fatal
	* sample/contig mismatch errors.
	*/
	void copy_header (const bcf_hdr_t *src) {
		// adds header lines from source header
		bcf_hdr_merge(hdr, src);
		// copy over samples manually
		for (int n = 0; n < bcf_hdr_nsamples(src); ++n) {
			bcf_hdr_add_sample(hdr, src->samples[n]);
		}
	}

	// Adds the tag from the given filter to the file header
	// This must be done before marking records with that filter.
	void add_filter_tag (const base_orient_bias_filter_f &filter) {
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