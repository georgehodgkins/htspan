#ifndef _HTSPAN_PRINT_HELP_HPP_
#define _HTSPAN_PRINT_HELP_HPP_

#include <stdexcept>
#include <iostream>
#include <cstring>

#include "optionparser.hpp"
#include "options.hpp"
#include "cstring.hpp"

/**
* This file contains some generic methods for printing subsets
* of the option list (not supported by native printUsage from the
* option parser).
*
* It also contains specific uses of those methods for hts-orient-bias.
*/

namespace hts {

namespace frontend {

using namespace option;

// Formatting markers to turn bold text on and off
const char* BOLD = "\e[1m";
const char* NO_BOLD = "\e[0m";

/**
* Print the help text of a given Descriptor,
* with word wrapping and indentation.
*
* @param descr The Descriptor to print
* @param prn The ostream to print to
* @param width Maximum number of columns written before wrapping.
* @return none (Prints to the selected ostream)
*/
void print_descr_usage(const Descriptor descr, std::ostream& prn, size_t width = 80) {
	size_t wctr = 0;
	bool tab_in = false;
	prn << BOLD;
	for (size_t n = 0; n < strlen(descr.help); ++n) {
		// CR indicates the end of the title
		if (descr.help[n] == '\r') {
			tab_in = true;
			wctr = 0;
			prn << NO_BOLD << "\n\t";
		} else if (descr.help[n] == '\n') {
			wctr = 0;
			prn << '\n';
			if (tab_in) {
				prn << '\t';
			}
		// whole-word wrapping: wraps if the next word (series of chars until a space or string end)
		// is longer than the remaining width in this row
		} else if (descr.help[n] == ' ') {// space, 0x20 
			const char *next_space = strchr(&descr.help[n+1], ' ');
			if ((next_space == NULL && strlen(&descr.help[n+1]) > (width - wctr)) ||
			(next_space > &descr.help[n + (width - wctr)])) {
				wctr = 0;
				prn << '\n';
				if (tab_in) {
					prn << '\t';
				}
			} else {
				prn << ' ';
			}
		} else {
			prn << descr.help[n];
			++wctr;
		}
	}
	prn << "\n\n";
}

/**
* This function prints a subset of help texts given
* the enum indices of the desired options.
*
* The complete usage array must be indexed with signed enums, 
* and must be sorted in ascending order. This is because multiple
* options can share the same enum index, so an option may not be at
* the position indicated by its index in the descriptor array,
* so it can be necessary to search for it (and those restrictions
* greatly simplify searching).
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
		// search for the correct option if it is not right at its index in the usage array
		// only works if usage array is sorted
		if (full_usage[i].index > selection[x]) {
			while (full_usage[i].index != selection[x] && i >= 0) {
				--i;
			}
		} else if (full_usage[i].index < selection[x]) {
			while (full_usage[i].index != selection[x] && i < (signed int) n_full) {
				++i;
			}
		}
		// if i goes below zero, it will wrap around and be > n_full
		if (i >= (signed int) n_full) {
			throw std::runtime_error("Error: a requested option was not found in the usage list (check that the list is sorted).");
		}
		print_descr_usage(full_usage[i], prn);
	}
}

/**
* Header printed at the top of all help texts
*/
void print_help_info () {
	std::cerr << BOLD << 
	"\nhtspan orient-bias: " << NO_BOLD << "Corrects errors in high-throughput sequencing data.\n" <<
	"Values in [brackets] are defaults.\n\n";
}

/**
* Print a list of available commands.
*/
void print_general_usage() {
	print_help_info();
	std::cerr << 
		"Commands:\n" <<
		BOLD << "quantify " << NO_BOLD << "- Calculate an estimate of global damage in the given BAM file.\n" <<
		BOLD << "identify " << NO_BOLD << "- Examine each SNV in the given snv file and provide a likelihood of damage.\n\n" <<

		"Use --help [command]  or --help [option] (no dashes) for details on use\n" <<
		"or --help utility for information on options relating to logging, debugging, and output.\n";
	std::cerr.flush();
}

// The various arrays which specify which options to print are in options.hpp

/**
* Print the usage information for the 'quantify' command.
*/
void print_quant_usage() {
	print_help_info();
	std::cerr << BOLD << "Quantification help:\n\n" << NO_BOLD;
	print_selected_usages (usage, total_arg_count, common_args, common_arg_count);
	print_selected_usages (usage, total_arg_count, quant_only_args, quant_arg_count);
}

/**
* Print the usage information for the 'identify' command.
*/
void print_ident_usage() {
	print_help_info();
	std::cerr << BOLD << "\nIdentification help:\n" << NO_BOLD;
	print_selected_usages (usage, total_arg_count, common_args, common_arg_count);
	print_selected_usages (usage, total_arg_count, ident_only_args, ident_arg_count);
}

/**
* Print the usage information for misc utility commands.
*/
void print_utility_usage() {
	print_help_info();
	std::cerr << BOLD << "Utilty flags:\n" << NO_BOLD;
	print_selected_usages (usage, total_arg_count, utility_args, utility_arg_count);
}

/*
* Print the correct help message based on the argument passed to --help.
*/
void print_help (const char* arg) {
	if (arg == NULL) {
		print_general_usage();
	} else if (strcmpi(arg, "quantify") == 0) {
		print_quant_usage();
	} else if (strcmpi(arg, "identify") == 0) {
		print_ident_usage();
	} else if (strcmpi(arg, "utility") == 0) {
		print_utility_usage();
	} else {
		if (strlen(arg) == 1) {
			for (size_t n = 0; n < total_arg_count; ++n) {
				if (strspn(arg, usage[n].shortopt) > 0) {
					print_descr_usage(usage[n], std::cerr);
					return;
				}
			}
		} else {
			for (size_t n = 0; n < total_arg_count; ++n) {
				if (strcmpi(arg, usage[n].longopt) == 0) {
					print_descr_usage(usage[n], std::cerr);
					return;
				}
			}
		}
		if (strlen(arg) > 0) {
			std::cerr << "Unknown help option \'" << arg << "\'.\n\n";
		}
		print_general_usage();
	}
}

}// namespace frontend

}// namespace hts

#endif // _HTSPAN_PRINT_HELP_HPP_