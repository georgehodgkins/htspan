#ifndef _HTSPAN_SNV_WRITER_HPP_
#define _HTSPAN_SNV_WRITER_HPP_

#include <stdexcept>
#include <cstring>

#include "htslib/vcf.h"

namespace hts {

namespace snv {

using namespace std;

// TODO: documentation

struct vcf_writer {
	// Pointer to VCF file
	htsFile *hf;
	// Pointer to BCF header struct for output file
	bcf_hdr_t *hdr;

	vcf_writer(const char* path, htsFile *f_templ, bcf_hdr_t *h_templ) {
		open(path, f_templ, h_templ);
	}

	void open(const char* path, htsFile *f_templ, bcf_hdr_t *h_templ) {
		hdr = bcf_hdr_dup(h_templ);
		htsFormat form = f_templ->format;
		hf = hts_open_format(path, "w", &form);
		if (!hf) {
			throw runtime_error("Error opening VCF file for output.");
		}
		int status = bcf_hdr_write(hf, hdr);
		if (status != 0) {
			throw runtime_error("Error writing header to VCF file.");
		}
	}

	void write(bcf1_t* rec, const bcf_hdr_t *src) {
		bcf_hdr_merge(hdr, src);
		int status = bcf_write(hf, hdr, rec);
		if (status != 0) {
			throw runtime_error("Error writing record to VCF file.");
		}
	}

	void close() {
		hts_close(hf);
		bcf_hdr_destroy(hdr);
		hdr = NULL;
		hf = NULL;
	}

	~vcf_writer () {
		close();
	}
};

}// namespace snv

}// namespace hts

#endif // _HTSPAN_SNV_WRITER_HPP_