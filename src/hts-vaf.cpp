#include <iostream>
#include <cstdlib>

#include "htspan/fetcher.hpp"
#include "htspan/nucleotide.hpp"
#include "htspan/orient_bias_data.hpp"
#include "htspan/io/snv_reader.hpp"
#include "htspan/io/snv_writer.hpp"

#include "frontend/orient-bias-identify.hpp"

// Quick-and-dirty implementation of a naive variant frequency filter

using namespace std;

int main (int argc, char** argv) {
	if (argc < 6) {
		cerr << "Usage: " << argv[0] << " <ref> <alt> <bam_file> <snv_file> <out_file> [<verbose>]";
	}

	nuc_t ref = char_to_nuc(argv[1][0]);
	nuc_t alt = char_to_nuc(argv[2][0]);
	string bam_file = argv[3];
	string snv_file = argv[4];
	string out_file = argv[5];

	hts::snv::tsv_reader snvr (snv_file.c_str());
	hts::snv::record rec;

	hts::snv::tsv_writer snvw (out_file.c_str(), hts::snv::F_TSV);
	snvw.add_numeric_info("VAFF", "Frequency of variant at the given locus");

	hts::fetcher f;
	if (!f.open(bam_file.c_str())) {
		cerr << "Could not open BAM file " << bam_file << endl;
		return 1;
	}

	hts::orient_bias_data data (ref, alt, 0); 	
	//cout << "snv\tfreq\n";
	while (fetch_next_snv(snvr, f, data, rec)) {
		if (snvr.error() == -1) {
			snvw.write(rec);
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
		//cout << rec.to_string() << '\t' << freq << '\n';
	}
	return 0;
}