#include <iostream>
#include <cstdlib>

#include "htspan/fetcher.hpp"
#include "htspan/nucleotide.hpp"
#include "htspan/orient_bias_data.hpp"
#include "htspan/io/snv_reader.hpp"
#include "htspan/io/snv_writer.hpp"

#include "frontend/orient-bias-identify.hpp"

/**
* Basic implementation of a naive variant frequency counter
* Counts SNV frequency for a specific variant type (and its complement)
* at each locus and records it in the output file.
*
* @arg ref Reference nuc for the variant
* @arg alt Alternate nuc for the variant
* @arg bam_file BAM file containing alignment of interest
* @arg snv_file VCF or TSV file containing variants to analyze (type autodetected)
* @arg out_file VCF or TSV file to write results to (type autodetected, no TSV->VCF conversion)
* @arg verbose Whether to print a table to stdout
*/

using namespace std;

int main (int argc, char** argv) {
	if (argc < 6) {
		cerr << "Usage: " << argv[0] << " <ref> <alt> <bam_file> <snv_file> <out_file> [<verbose>]\n";
		return 1;
	}

	nuc_t ref = char_to_nuc(argv[1][0]);
	nuc_t alt = char_to_nuc(argv[2][0]);
	string bam_file = argv[3];
	char* snv_file = argv[4];
	char* out_file = argv[5];
	bool verbose = false;
	if (argc > 6) {
		verbose = atoi(argv[6]) != 0;
	}

	hts::snv::streamer s (snv_file, out_file, hts::snv::opts_to_fmt(snv_file, NULL), hts::snv::opts_to_fmt(out_file, NULL));

	hts::snv::reader &snvr = *s.snvr_pt;
	hts::snv::writer &snvw = *s.snvw_pt;

	snvw.add_numeric_info("VAFF", "Frequency of variant at the given locus");

	hts::snv::record rec;

	hts::fetcher f;
	if (!f.open(bam_file.c_str())) {
		cerr << "Could not open BAM file " << bam_file << endl;
		return 1;
	}

	hts::orient_bias_data data (ref, alt, 0); 	
	
	if (verbose) cout << "snv\tfreq\n";

	while (fetch_next_snv(snvr, f, data, rec)) {
		if (snvr.error()) {
			if (snvr.error() == -1) { // successful read, inconsistent variant
				snvw.write(rec);
			}
			continue;
		}

		size_t count = 0;
		for (size_t i = 0; i < data.bases.size(); ++i) {
			if (data.bases[i] == 1) {
				++count;
			}
		}

		double freq = (double) count / data.bases.size();

		snvw.write(rec, "VAFF", freq);
		if (verbose) cout << rec.to_string() << '\t' << freq << '\n';

	}
	return 0;
}