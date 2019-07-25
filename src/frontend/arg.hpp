#ifndef _HTSPAN_ARG_HPP_
#define _HTSPAN_ARG_HPP_

#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <climits>
#include <cfloat>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "optionparser.hpp"
#include "../htspan/nucleotide.hpp"
#include "../htspan/piler.hpp"
#include "file.hpp"
#include "cstring.hpp"

// This file contains argument checks for the options defined in options.hpp

namespace hts {
	
namespace frontend {

using namespace option;

struct Arg: public option::Arg {
	//
	// Generic option checks used for multiple options
	//

	/**
	* Check whether file can be opened for input.
	* 
	* @param opt Option structure containing the potential filename
	* @param msg Whether to print a message upon rejection
	* @param xtns Extensions to require, in all lowercase (case insensitive)
	* @paran n_xtns Number of members in xtns
	* @return ARG_OK if filename passes and ARG_ILLEGAL if it does not
	*/
	static ArgStatus InputFile (const Option& opt, bool msg, const char* xtns[] = NULL, const size_t n_xtns = 0) {
		if (opt.arg == NULL || strlen(opt.arg) == 0) {
			if (msg) {
				std::cerr << "Option " << opt.name << " requires an argument. Try --help.\n";
			}
			return ARG_ILLEGAL;
		}
		if (!file_exists(opt.arg)) {
			if (msg) {
				std::cerr << "Could not open path \'" << opt.arg << "\' for reading (passed to " << opt.name << ").\n";
			}
			return ARG_ILLEGAL;
		}
		if (xtns != NULL && n_xtns > 0) {
			// find the beginning of the last extension in filename
			const char* dot = strrchr(opt.arg, '.') + 1;
			for (size_t i = 0; i < n_xtns; ++i) {
				if (strcmpi(dot, xtns[i]) == 0) {
					return ARG_OK;
				}
			}
			//if this point is reached no extension was matched
			if (msg) {
				std::cerr << "File passed to " << opt.name << " must have one of the following extensions: ";
				for (size_t n = 0; n < n_xtns-1; ++n) {
					std::cerr << '.' << xtns[n] << ", ";
				}
				std::cerr << '.' << xtns[n_xtns-1] << ".\n";
			}
			return ARG_ILLEGAL;
		} else { // n_xtns == 0, no extension checking
			return ARG_OK;
		}
	}

	/**
	* Check whether file can be opened for output.
	*
	* Params are the same as above, except with no extension checking
	*/
	static ArgStatus OutputFile (const Option& opt, bool msg) {
		if (opt.arg == NULL || strlen(opt.arg) == 0) {
			if (msg) {
				std::cerr << "Option " << opt.name << " requires an argument. Try --help.\n";
			}
			return ARG_ILLEGAL;
		}
		if (strlen(opt.arg) == 1 && opt.arg[0] == '-') {
			return ARG_OK;
		}
		if (!file_writable(opt.arg)) {
			if (msg) {
				std::cerr << "Could not open path \'" << opt.arg << "\' for writing (passed to " << opt.name << ").\n";
			}
			return ARG_ILLEGAL;
		} else {
			return ARG_OK;
		}
	}

	/**
	* Check that a nucleotide argument is canonical.
	*/
	static ArgStatus CanonicalNucleotide (const Option& opt, bool msg) {
		if (opt.arg == NULL || strlen(opt.arg) == 0) {
			if (msg) {
				std::cerr << "Option " << opt.name << " requires an argument. Try --help.\n";
			}
			return ARG_ILLEGAL;
		}
		ArgStatus rtn;
		if (strlen(opt.arg) != 1) {
			rtn = ARG_ILLEGAL;
		} else {
			rtn = (nuc_is_canonical(char_to_nuc(opt.arg[0]))) ? ARG_OK : ARG_ILLEGAL;
		}
		if (msg && rtn == ARG_ILLEGAL) {
			std::cerr << '\'' << opt.arg << "\' does not specify a canonical nucleotide (A, C, G, or T).\n";
		}
		return rtn;
	}

	/**
	* Check that a integral argument is valid and in the range [lb, ub].
	*/
	static ArgStatus IntRange (const Option& opt, bool msg, int lb = INT_MIN, int ub = INT_MAX) {
		if (opt.arg == NULL || strlen(opt.arg) == 0) {
			if (msg) {
				std::cerr << "Option " << opt.name << " requires an argument. Try --help.\n";
			}
			return ARG_ILLEGAL;
		}
		int x;
		// strtol returns zero for an invalid argument, so inputs of zero require special handling
		// check if the arg string is all zeroes and convert without strtol if so
		if (strspn(opt.arg, "0") == strlen(opt.arg)) {
			x = 0;
		} else {
			x = strtol(opt.arg, NULL, 10);
			// here, x == 0 indicates the input does not describe a valid integer
			if (x == 0) {
				if (msg) {
					std::cerr << "\'" << opt.arg << "\' does not specify a valid integer (passed to " << opt.name << ").\n";
				}
				return ARG_ILLEGAL;
			}
		}
		
		// bounds checking
		if (x <= ub && x >= lb) {
			return ARG_OK;
		} else {
			if (msg) {
				std::cerr << "Argument to " << opt.name << " must fall between " << lb << " and " << ub << ", inclusive.\n";
			}
			return ARG_ILLEGAL;
		}
	}

	/**
	* Check that a double argument is valid and in the range [lb, ub]
	*/
	static ArgStatus DoubleRange (const Option& opt, bool msg, double lb = DBL_MIN, double ub = DBL_MAX) {
		if (opt.arg == NULL || strlen(opt.arg) == 0) {
			if (msg) {
				std::cerr << "Option " << opt.name << " requires an argument. Try --help.\n";
			}
			return ARG_ILLEGAL;
		}
		// see int check function above for explanation of the zero-handling
		double x;
		// any valid, normal way to input a zero fp number should match a substring of this string
		// corner case: this also matches ".", which is not a valid input, hence the second condition
		if (strstr("000.000", opt.arg) != NULL && (opt.arg[0] != '.' || strlen(opt.arg) > 1)) {
			x = 0.0;
		} else {
			x = strtod(opt.arg, NULL);
			if (x == 0.0) {
				if (msg) {
					std::cerr << "\'" << opt.arg << "\' does not specify a valid floating-point value (passed to " << opt.name << ").\n";
				}
				return ARG_ILLEGAL;
			}
		}
		
		// bounds checking
		if (x <= ub && x >= lb) {
			return ARG_OK;
		} else {
			if (msg) {
				std::cerr << "Argument to " << opt.name << " must fall between " << lb << " and " << ub << ", inclusive.\n";
			}
			return ARG_ILLEGAL;
		}
	}

	/**
	* Check that a string argument matches one of a list of options.
	*/
	static ArgStatus MatchingString (const Option& opt, bool msg, const char* str[], const size_t n_str) {
		for (size_t n = 0; n < n_str; ++n) {
			if (strcmpi(opt.arg, str[n]) == 0) {
				return ARG_OK;
			}
		}
		
		// if we reach this point arg did not match any valid strings
		if (msg) {
			std::cerr << "Argument to " << opt.name << " must match one of the following:\n";
			for (size_t n = 0; n < n_str-1; ++n) {
				std::cerr << "\'" << str[n] << "\', ";
			}
			std::cerr << "\'" << str[n_str-1] << "\'.\n";
		}
		return ARG_ILLEGAL;
	}

	//
	// Option-specific calls to generic checks
	//

	static ArgStatus DamageType (const Option& opt, bool msg) {
		const char* str[] = {"ffpe", "oxog"};
		return MatchingString(opt, msg, str, sizeof(str)/sizeof(char*));
	} 

	static ArgStatus ExternalSim (const Option& opt, bool msg) {
		const char* xtns[] = {"tsv"};
		return InputFile(opt, msg, xtns, sizeof(xtns)/sizeof(char*));
	}

	static ArgStatus Verbosity (const Option& opt, bool msg) {
		return IntRange (opt, msg, 0, 3);
	}


	static ArgStatus AlignmentFile (const Option& opt, bool msg, bool try_open = true) {
		// normal file and extension checks
		const char* xtns[] = {"bam"};
		ArgStatus rtn = InputFile(opt, msg, xtns, sizeof(xtns)/sizeof(char*));
		if (rtn != ARG_OK) {
			return rtn;
		}

		// check if the file is a valid BAM file with header and index
		if (try_open) {
			hts::piler p;
			if (!p.open(opt.arg)) {
				if (msg) {
					std::cerr << "Could not open \'" << opt.arg << "\' as an alignment file.";
				}
				return ARG_ILLEGAL;
			}
		}
		return ARG_OK;
	}

	static ArgStatus PairedEndAlignmentFile (const Option& opt, bool msg) {
		// number of reads to check
		const int n_checks = 10;

		// this function will open the file in the piler if it can be opened
		ArgStatus rtn = AlignmentFile (opt, msg, false);
		if (rtn != ARG_OK) {
			return rtn;
		}

		// open the file
		hts::piler p;
		if (!p.open(opt.arg)) {
			if (msg) {
				std::cerr << "Could not open \'" << opt.arg << "\' as an alignment file.";
			}
			return ARG_ILLEGAL;
		}

		// check up to 10 reads
		int total = 0;
		while (total < n_checks) {
			// get some reads
			const vector<bam1_t*> &pile = p.next();
			if (pile.empty()) continue; // no reads passing the qfilter, keep going
			if (pile[0] == NULL) break; // error or EOF
			
			// if any read in the pile has a read1/read2 flag, consider it paired-end
			for (size_t r = 0; r < pile.size(); ++r) {
				if (bam_is_read1 (pile[r]) || bam_is_read2 (pile[r])) {
					return ARG_OK;
				}
				++total;
			}
		}

		// none of the checked reads are paired-end: error
		if (msg) {
			std::cerr << "BAM file \'" << opt.arg << "\' is not paired-end; "
				"it appears to be single-end.";
		}
		return ARG_ILLEGAL;
	}

	static ArgStatus Model (const Option& opt, bool msg) {
		const char* str[] = {"bayes", "freq"};
		return MatchingString(opt, msg, str, sizeof(str)/sizeof(char*));
	}

	static ArgStatus ReferenceFile (const Option& opt, bool msg) {
		const char* xtns[] = {"fasta", "fa"};
		return InputFile(opt, msg, xtns, sizeof(xtns)/sizeof(char*));
	}

	static ArgStatus InSnvFile (const Option& opt, bool msg) {
		const char* xtns[] = {"tsv", "snv", "vcf", "gz", "bgz"};
		return InputFile(opt, msg, xtns, sizeof(xtns)/sizeof(char*));
	}

	static ArgStatus Probability (const Option& opt, bool msg) {
		return DoubleRange(opt, msg, 0, 1);
	}

	static ArgStatus MinMapQ (const Option& opt, bool msg) {
		return IntRange(opt, msg, 0, 100);
	}

	static ArgStatus MinBaseQ (const Option& opt, bool msg) {
		return IntRange(opt, msg, 0, 100);
	}

	static ArgStatus MaxQReads (const Option& opt, bool msg) {
		return IntRange(opt, msg, 100, INT_MAX);
	}

	static ArgStatus MinzBound (const Option& opt, bool msg) {
		return IntRange(opt, msg, 5, 1000);
	}

	static ArgStatus Eps (const Option& opt, bool msg) {
		return DoubleRange(opt, msg, 1e-12, .01);
	}

	static ArgStatus MinzIter (const Option& opt, bool msg) {
		return IntRange(opt, msg, 10, 1000);
	}

	static ArgStatus PositiveDouble (const Option& opt, bool msg) {
		return DoubleRange(opt, msg, 0.0, DBL_MAX);
	}

	static ArgStatus SnvFType (const Option& opt, bool msg) {
		const char* str[] = {"tsv", "vcf", "vcf-bgz"};
		return MatchingString(opt, msg, str, sizeof(str)/sizeof(char*));
	}

};// struct Arg

}// namespace frontend

}// namespace hts

#endif // _HTSPAN_ARG_HPP_
