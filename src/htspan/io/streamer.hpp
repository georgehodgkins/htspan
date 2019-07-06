#ifndef _HTSPAN_STREAMER_HPP_
#define _HTSPAN_STREAMER_HPP_

#include "snv.hpp"
#include "snv_reader.hpp"
#include "snv_writer.hpp"

namespace hts {

namespace snv {

struct streamer {
	reader *snvr_pt;
	writer *snvw_pt;

	/**
	* This function constructs the necessary SNV reader and writer
	* given a format and paths for the files, and returns base class
	* references to the constructed objects.
	*
	*/
	streamer(const char* read_path, const char* write_path, FMTFLAGS_T read_fmt, FMTFLAGS_T write_fmt) {
		using namespace snv;
		if ((read_fmt & F_TSV) && (write_fmt & F_VCF)) {
			throw runtime_error("Converting from TSV to VCF is not supported.");
		}
		// setup reader
		if (read_fmt & F_TSV) {
			snvr_pt = (reader*) new tsv_reader(read_path);
		} else if (read_fmt & F_VCF) {
			snvr_pt = (reader*) new vcf_reader(read_path);
		}
		// setup writer
		if (write_fmt & F_TSV) {
			snvw_pt = (writer*) new tsv_writer (write_path, write_fmt);
		} else if (write_fmt & F_VCF) {
			vcf_writer *vcfw_pt = new vcf_writer (write_path, write_fmt);
			vcfw_pt->copy_header(((vcf_reader*) snvr_pt)->hdr);
			snvw_pt = (writer*) vcfw_pt;
		}
	}

	~streamer() {
		delete snvr_pt;
		delete snvw_pt;
	}
};

} // namespace snv

} // namespace hts

#endif // _HTSPAN_STREAMER_HPP_