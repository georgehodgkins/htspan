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

#ifndef _HTSPAN_OPTCHK_HPP_
#define _HTSPAN_OPTCHK_HPP_

// UNDER CONSTRUCTION
// This file contains argument checks for the option parser.

namespace hts {
	
namespace frontend {


/**
* Case insensitive strcmp, returns bool instead of the traditional int
*/
int strcmpi (const char* s1, const char* s2) {
	while (tolower(*s1) == tolower(*s2)) {
		++s1;
		++s2;
		if (*s1 == 0x0) break;
	}
	if (*s1 == *s2) {//they should both be null so no need for tolower
		return 0;
	} else {
		return *s1 - *s2;
	}
}

using namespace option;

struct Arg: public option::Arg {
	/**
	* Check whether file can be opened for input.
	* 
	* @param opt Option structure containing the potential filename
	* @param msg Whether to print a message upon rejection
	* @param xtns Extensions to require, in all lowercase (case insensitive)
	* @paran n_xtns Number of members in xtns
	* @return ARG_OK if filename passes and ARG_ILLEGAL if it does not
	*/
	static ArgStatus InputFile (const Option& opt, bool msg, const char* xtns[], const size_t n_xtns = 0) {
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
		if (n_xtns > 0) {
			const char* dot = strrchr(opt.arg, '.') + 1;
			for(size_t i = 0; i < n_xtns; ++i) {
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
	* Check that a nucleotide argument is canonical
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
	* Check that a integral argument is valid and in the range [lb, ub]
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
	*
	* Since most of the arguments this will handle are probabilities,
	* lb and ub default to 0 and 1 respectively.
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
		if (x <= ub && x >= lb) {
			return ARG_OK;
		} else {
			if (msg) {
				std::cerr << "Argument to " << opt.name << " must fall between " << lb << " and " << ub << ", inclusive.\n";
			}
			return ARG_ILLEGAL;
		}
	}

	static ArgStatus DamageType (const Option& opt, bool msg) {
		if (strcmpi(opt.arg, "ffpe") == 0 ||
			strcmpi(opt.arg, "oxog") == 0) {
			return ARG_OK;
		} else {
			if (msg) {
				std::cerr << "Argument to " << opt.name << " must be either \'ffpe\' or \'oxog\'.";
			}
			return ARG_ILLEGAL;
		}
	} 

	static ArgStatus ExternalSim (const Option& opt, bool msg) {
		const char* xtns[] = {"tsv"};
		return InputFile(opt, msg, xtns, 1);
	}

	static ArgStatus Verbosity (const Option& opt, bool msg) {
		return IntRange (opt, msg, 0, 3);
	}


	static ArgStatus AlignmentFile (const Option& opt, bool msg, hts::piler& p) {
		// normal file and extension checks
		const char* xtns[] = {"bam"};
		ArgStatus rtn = InputFile(opt, msg, xtns, 1);
		if (rtn != ARG_OK) {
			return rtn;
		}
		// check if the file is a valid BAM file with header and index
		if (!p.open(opt.arg)) {
			if (msg) {
				std::cerr << "Could not open \'" << opt.arg << "\' as an alignment file.";
			}
			return ARG_ILLEGAL;
		}
		return ARG_OK;
	}

	static ArgStatus PairedEndAlignmentFile (const Option& opt, bool msg) {
		hts::piler p;
		// this function will open the file in the piler if it can be opened
		ArgStatus rtn = AlignmentFile (opt, msg, p);
		if (rtn != ARG_OK) {
			return rtn;
		}
		// read the first pile from the file
		const bam_pileup1_t *samp = p.next();
		if (samp == NULL) {
			if (msg) {
				std::cerr << "Could not get a read from \'" << opt.arg << "\'.";
			}
			return ARG_ILLEGAL;
		}
		// if any read in the pile has a read1/read2 flag, consider it paired-end
		// otherwise, error (after 10 tests)
		int rfs = 0; // counts reads from the current pileup
		for (size_t r = 0; r < 10; ++r) {
			// handles cases where the pileup has <10 reads
			if (rfs >= p.curr_plp_size()) {
				samp = p.next();
				rfs = 0;
			}
			if (bam_is_read1 (samp[r].b) || bam_is_read2 (samp[r].b)) {
				p.close();
				return ARG_OK;
			}
			++rfs;
		}
		p.close();
		if (msg) {
			std::cerr << "BAM file \'" << opt.arg << "\' is not paired end.";
		}
		return ARG_ILLEGAL;
	}


	static ArgStatus ReferenceFile (const Option& opt, bool msg) {
		const char* xtns[] = {"fasta", "fa"};
		return InputFile(opt, msg, xtns, 2);
	}

	static ArgStatus SnvFile (const Option& opt, bool msg) {
		const char* xtns[] = {"tsv", "snv"};
		return InputFile(opt, msg, xtns, 2);
	}

	static ArgStatus Probability (const Option& opt, bool msg) {
		return DoubleRange(opt, msg, 0, 1);
	}

	static ArgStatus MinMapQ (const Option& opt, bool msg) {
		return IntRange(opt, msg, 0, 100);// check bounds
	}

	static ArgStatus MinBaseQ (const Option& opt, bool msg) {
		return IntRange(opt, msg, 0, 100);// ^^^
	}

	static ArgStatus MaxQReads (const Option& opt, bool msg) {
		return IntRange(opt, msg, 1e6, 1e9);
	}

	static ArgStatus MinzBound (const Option& opt, bool msg) {
		return IntRange(opt, msg, 5, 1000);
	}

	static ArgStatus MinzEps (const Option& opt, bool msg) {
		return DoubleRange(opt, msg, 1e-12, .01);
	}

	static ArgStatus MinzIter (const Option& opt, bool msg) {
		return IntRange(opt, msg, 10, 1000);
	}

	static ArgStatus ErrMeanSim (const Option& opt, bool msg) {// check bounds
		return DoubleRange(opt, msg, 0.0, 100.0);
	}

	static ArgStatus ErrStdevSim (const Option& opt, bool msg) {// ^^^
		return DoubleRange(opt, msg, 0.0, 100.0);
	}

};// struct Arg

}// namespace frontend

}// namespace hts

#endif // _HTSPAN_OPTCHK_HPP_