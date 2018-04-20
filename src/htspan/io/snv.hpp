#ifndef _HTSPAN_SNV_HPP_
#define _HTSPAN_SNV_HPP_

#include <fstream>
#include <sstream>
#include <stdexcept>

#include "htspan/nucleotide.hpp"


namespace hts {

using namespace std;

namespace snv {

struct record {
	string chrom;
	int32_t pos;
	nuc_t nt_ref;
	nuc_t nt_alt;
};

struct reader {
	ifstream f;

	string line;

	reader(const char* path)
	: f(path) {
		if (!f.is_open()) {
			throw runtime_error("Error: could not open input file");
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

		// convert from 1-based to 0-based
		r.pos -= 1;
		r.nt_ref = char_to_nuc(char_ref);
		r.nt_alt = char_to_nuc(char_alt);

		return true;
	}

	void close() {
		f.close();
	}
};

}  // namespace snv

}  // namepsace hts

#endif  // _HTSPAN_SNV_HPP_
