#include <stdexcept>
#include <iostream>
#include <cstring>

#include "optionparser.hpp"
#include "options.hpp"

namespace hts {

namespace frontend {

using namespace option;

/**
* Print the help text of a given descriptor.
* Needed because the one that comes with the parser
* requires everything to be const.
*
* @param descr The Descriptor to print
* @param prn The ostream to print to
* @param width Maximum number of columns written before wrapping.
* @return none (Prints to the selected ostream)
*/
void print_descr_usage(const Descriptor descr, std::ostream& prn, int width = 80) {
	int wctr = 0;
	for (size_t n = 0; n < strlen(descr.help); ++n) {
		if (wctr == width) {
			prn << std::endl;
			wctr = 0;
		} else if (descr.help[n] == '\n') {
			wctr = -1;
		}
		prn << descr.help[n];
		++wctr;
	}
	prn << std::endl;
}

/**
* This function prints a subset of help texts given
* the enum indices of the desired options.
*
* @param full_usage The complete usage array
* @param n_full Length of full_usage
* @param selection Array of enum indices whose corresponding help texts should be printed
* @param n_sel Length of selection
* @param prn The ostream to write to [cerr]
* @return none (Writes to some ostream, probably cerr)
*/
void print_selected_usages(const Descriptor full_usage[], size_t n_full, const OptionIndex selection[], const size_t n_sel, std::ostream& prn = std::cerr) {
	for (size_t x = 0; x < n_sel; ++x) {
		signed int i = selection[x];//must be signed so that lower bound checking works
		// search for the correct option if it is not right at its index
		// only works if usage array is sorted
		if (full_usage[selection[x]].index > selection[x]) {
			while (full_usage[i].index != selection[x] && i >= 0) {
				--i;
			}
		} else if (full_usage[selection[x]].index < selection[x]) {
			while (full_usage[i].index != selection[x] && i < n_full) {
				++i;
			}
		}
		if (i < 0 || i >= n_full) {
			throw std::runtime_error("Error: a requested option was not found in the usage list (check that the list is sorted).");
		}
		print_descr_usage(full_usage[i], prn);
	}
	prn << std::endl; // print an extra line after the full usage block
}

void print_general_usage() {
	std::cerr << 
		"htspan: Corrects errors in high-throughput sequencing data.\n" <<
		"Help format: -f --flag usage [default]\n" <<
		"more usage\n\n" <<

		"Commands:\n" <<
		"quantify - Calculate an estimate of global damage in the given BAM file.\n" <<
		"identify - Examine each SNV in the given snv file and provide a likelihood of damage.\n\n" <<

		"Use --help [command] for details on use or\n" <<
		"--help utility for information on options relating to logging, debugging, and output.\n";
	std::cerr.flush();
}

// The various arrays which specify which options to print are in options.hpp

/**
* Print the usage information for the 'quantify' command.
*/
void print_quant_usage() {
	std::cerr << "\nQuantification help:\n" <<
		"Common flags:\n";
	print_selected_usages (usage, total_arg_count, common_args, common_arg_count);
	std::cerr << "Quantification-specific flags:\n";
	print_selected_usages (usage, total_arg_count, quant_only_args, quant_arg_count);
}

/**
* Print the usage information for the 'identify' command.
*/
void print_ident_usage() {
	std::cerr << "\nIdentification help:\n" << 
		"Common flags:\n";
	print_selected_usages (usage, total_arg_count, common_args, common_arg_count);
	std::cerr << "Identification-specific flags:\n";
	print_selected_usages (usage, total_arg_count, quant_only_args, quant_arg_count);
}

/**
* Print the usage information for misc utility commands.
*/
void print_utility_usage() {
	std::cerr << "Utilty flags:\n";
	print_selected_usages (usage, total_arg_count, utility_args, utility_arg_count);
}

/*
* Print the correct help message based on the argument given.
*/
void print_help (const char* arg) {
	if (arg == NULL) {
		print_general_usage();
	} else if (strcmp(arg, "quantify") == 0) {
		print_quant_usage();
	} else if (strcmp(arg, "identify") == 0) {
		print_ident_usage();
	} else if (strcmp(arg, "utility") == 0) {
		print_utility_usage();
	} else {
		if (strlen(arg) > 0) {
			std::cerr << "Unknown help option \'" << arg << "\'.\n\n";
		}
		print_general_usage();
	}
}

}// namespace frontend

}// namespace hts