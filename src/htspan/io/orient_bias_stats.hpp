#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>

namespace hts {

namespace orient_bias_stats {
	using namespace std;
	struct record {
		int base;
		double error;
		bool orient;
	};

	struct reader {
		ifstream readin;
		size_t valid_lines;
		
		reader (const char* path) {
			string tmp;
			readin.open(path, ifstream::in);
			if (!readin.is_open()) {
				throw runtime_error("Error: could not open stat file");
			}
			//discard header line
			getline(readin, tmp);
		}

		//get next record
		bool next (record& rec) {
			string tmp;
			//ignore empty lines and comments marked with #
			while (readin.peek() == '\n' || readin.peek() == '#') {
				getline(readin, tmp);
			}
			//return false if done
			if (readin.eof()) {
				return false;
			}
			getline(readin, tmp);
			stringstream line (tmp);
			char r_ch;
			//get base (0=ref or 1=alt), error, and first char of orient bool
			line >> rec.base >> rec.error >> r_ch;
			rec.orient = (r_ch == 'T');//bools recorded as TRUE or FALSE
			//we read a line
			return true;
		}

		void close() {
			readin.close();
		}

	};
}//namespace orient_bias_stats

}//namespace hts