#ifndef _HTSPAN_SNV_WRITER_HPP_
#define _HTSPAN_SNV_WRITER_HPP_

#include <stdexcept>
#include <cstring>
#include <fstream>

#include "htslib/vcf.h"
#include "htslib/hts.h"

#include "snv_reader.hpp"
#include "../base_orient_bias_filter.hpp"

namespace hts {

namespace snv {

using namespace std;

// TODO: documentation

struct tsv_writer {
	ofstream f;

	tsv_writer(const char* path, tsv_reader &to_copy) {
		open(path);
	}

	void open (const char* path) {
		f.open(path, ios::trunc);
		if (!f.is_open()) {
			throw runtime_error("Error opening TSV file for output.");
		}
		// write header line
		f << "chrom\tpos\tref\talt" << endl;
	}

	void write (record r, const char* line = NULL) {
		f << r.chrom << '\t' << r.pos+1 << '\t' << nuc_to_char(r.nt_ref) << '\t' << nuc_to_char(r.nt_alt) << endl;
	}

	void close() {
		f.close();
	}

	void add_filter_tag (base_orient_bias_filter_f& filter) {}

	~tsv_writer() {
		close();
	}
};



struct vcf_writer {
	// Pointer to VCF file
	htsFile *hf;
	// Pointer to BCF header struct for output file
	bcf_hdr_t *hdr;
	// Linked header (to synchronize with
	bcf_hdr_t *linked;

	vcf_writer(const char* path, vcf_reader &to_copy) {
		open(path, to_copy.hf, to_copy.hdr);
	}

	void open(const char* path, htsFile *f_templ, bcf_hdr_t *h_templ) {
		linked = h_templ;
		hdr = bcf_hdr_dup(linked);
		char mode[2];
		mode [0] = 'w';
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
		int status = bcf_hdr_write(hf, hdr);
		if (status != 0) {
			throw runtime_error("Error writing header to VCF file.");
		}
	}

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
		int status = bcf_write(hf, hdr, rec);
		if (status != 0) {
			throw runtime_error("Error writing record to VCF file.");
		}
	}

	void sync_header() {
		bcf_hdr_merge(hdr, linked);
	}

	void add_filter_tag (base_orient_bias_filter_f& filter) {
		sync_header();
		bcf_hdr_printf(hdr, "FILTER=<ID=%s,Description=\"%s\">", filter.text_id, filter.get_description());
	}

	void close() {
		linked = NULL;
		if (hf != NULL) {
			hts_close(hf);
			hf = NULL;
		}
		if (hdr != NULL) {
			bcf_hdr_destroy(hdr);
			hdr = NULL;
		}
	}

	bool is_annotated () const {
		return true;
	}

	~vcf_writer () {
		close();
	}
};

}// namespace snv

}// namespace hts

#endif // _HTSPAN_SNV_WRITER_HPP_