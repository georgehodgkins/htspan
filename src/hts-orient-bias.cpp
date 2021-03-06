#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <utility>

#include "frontend/options.hpp"
#include "frontend/print-help.hpp"
#include "frontend/orient-bias-identify.hpp"
#include "frontend/orient-bias-quantify.hpp"
#include "frontend/optionparser.hpp"
#include "frontend/string.hpp"
#include "frontend/simul_writer.hpp"
#include "frontend/json.h"

#include "htspan/nucleotide.hpp"
#include "htspan/fetcher.hpp"
#include "htspan/piler.hpp"
#include "htspan/print.hpp"

#include "htspan/io/faidx_reader.hpp"
#include "htspan/io/snv.hpp"

using namespace hts::frontend;

/**
* This is the entry point for hts-orient-bias. The code in this file
* is mainly concerned with parsing command-line options
* and initializing input/output files; it then passes these to driver functions
* located in frontend/orient-bias-identify.hpp and frontend/orient-bias-quantify.hpp.
*/
int main (int argc, char** argv) {
	// enum for model options
	enum ModelType {BAYES, FREQ};

	//
	// Parameter defaults
	// The value of these is noted as a comment where they are used to set parameters;
	// if you change a value here, change it there too
	//
	const bool default_use_stdout = true;
	const int default_verbosity = 1;
	const string default_snv_out_fname = "out.snv";
	const ModelType default_model = FREQ;
	// default siglevel depends on the model used, since they return different statistics
	const double default_siglevel_freq = .05;
	const double default_siglevel_bayes = .95;
	const double default_phi = .01;
	const double default_alpha = 1.0;
	const double default_beta = 1.0;
	const double default_prior_alt = .5;
	const int default_min_mapq = 30;
	const int default_min_baseq = 20;
	const int default_max_reads = 1000000;
	const int default_minz_bound = 15;
	const double default_eps = 1e-6;
	const bool default_keep_dup = false;


	//
	// Process input arguments
	// 
	if (argc <= 1) { // no arguments given
		print_help(NULL);
		return 1;
	}
	argv++; // ignore executable name
	argc--;
	// first argument should be the command:
	// 'identify' or 'quantify', or a help flag
	const char* command = argv[0];
	// get command
	bool quantifying = false;
	bool identifying = false;
	if (strcmpi(command, "quantify") == 0) {
		quantifying = true;
		// skip command
		argv++;
		argc--;
	} else if (strcmpi(command, "identify") == 0) {
		identifying = true;
		// skip command
		argv++;
		argc--;
	} else if (strcmpi(command, "--help") == 0 ||
			strcmpi(command, "-?") == 0 ||
			strcmpi(command, "help") == 0) {
		
		print_help(argv[1]);
		return 0;
	}
	// parse options according to definitions in options.hpp
	option::Stats stats(usage, argc, argv);
	option::Option options[stats.options_max], buffer[stats.buffer_max];
	option::Parser parse(usage, argc, argv, options, buffer);
	if (parse.error()) {
		std::cerr << "Error encountered in parsing.\n";
		print_help(NULL);
		return 1;
	}
	// print the requested help message and quit if the help flag is passed
	if (options[HELP]) {
		if (options[HELP].arg) {
			print_help(options[HELP].arg);
			return 0;
		} else {
			print_help(command);
			return 0;
		}
	} else if (!quantifying && !identifying) {
		std::cerr << "\'" << command << "\' is not a recognized command.\n\n";
		print_help(NULL);
		return 1;
	}
	// Convert parsed options to runtime vars
	// Note that all arg checking is handled by the parser,
	// so args here are assumed to be individually valid,
	// although some combinations are still invalid

	//
	// Argument parsing block
	// 
	bool success = true;

	// Set up logging output
	bool use_stdout = default_use_stdout; // true
	if (options[STDOUT]) {
		use_stdout = options[STDOUT].last()->type() == t_ON;
	}
	int verbosity = default_verbosity; // 1
	if (options[VERBOSITY]) {
		verbosity = atoi(options[VERBOSITY].arg);
	}
	bool log_to_file = (bool) options[LOGFILE];
	std::string log_fname;
	if (log_to_file) {
		log_fname = options[LOGFILE].arg;
	}
	std::string json_fname;
	if (options[JSON_OUT]) {
		json_fname = options[JSON_OUT].arg;
	}

	// plain output makes output machine-readable
	// Errors and warnings will still be printed, though
	bool plain_output = (bool) options[PLAIN];

	// Direct file and result output appropriately
	
	// Fatal errors will go to stderr regardless either by
	// cerr in the frontend or thrown exceptions in the backend
	if (log_to_file) {
		global_log.add_file(log_fname);
	}
	global_log.use_cout(use_stdout);
	
	// Each of the driver functions has a section in its header listing
	// what is output at each verbosity level
	global_log.set_verbosity(verbosity);

	// Ref/alt nucleotide options (usually set by damage type flag)
	nuc_t ref = nuc_N;
	nuc_t alt = nuc_N;
	string dtype;
	if (options[DAMAGE_TYPE]) {
		if (strcmpi(options[DAMAGE_TYPE].arg, "ffpe") == 0) {
			ref = nuc_C;
			alt = nuc_T;
			dtype = "ffpe";
		} else if (strcmpi(options[DAMAGE_TYPE].arg, "oxog") == 0) {
			ref = nuc_G;
			alt = nuc_T;
			dtype = "oxog";
		}
	} else if (options[REF] && options[ALT]) {
		ref = char_to_nuc(options[REF].arg[0]);
		alt = char_to_nuc(options[ALT].arg[0]);
		dtype = "other";
		if (nuc_equal (alt, ref)) {
			std::cerr <<
				"Error: The reference and alternative nucleotides cannot be the same.\n";
			success = false;
		}
	} else {
		std::cerr <<
			"A variant type must be specified either with the -t/--damage-type flag or using the -R and -A flags.\n";
			success = false;
	}

	// Alignment file options
	std::string align_fname;
	if (!options[BAMFILE]) {
		std::cerr << 
			"Error: alignment file argument (-b, --alignment-file) is required.\n";
		success = false;
	} else {
		align_fname = options[BAMFILE].arg;
	}

	// Reference file options
	std::string ref_fname;
	if (options[REFFILE]) {
		ref_fname = options[REFFILE].arg;
	}

	// Input SNV file options
	std::string snv_in_fname;
	snv::FMTFLAGS_T snv_in_fmt = snv::F_NULL;
	if (options[IN_SNVFILE]) {
		snv_in_fname = options[IN_SNVFILE].arg;
		// opts_to_fmt deduces format either from an explicitly set option or from the filename
		if (options[IN_SNVFTYPE]) {
			snv_in_fmt = hts::snv::opts_to_fmt(snv_in_fname.c_str(), options[IN_SNVFTYPE].arg);
		} else {
			snv_in_fmt = hts::snv::opts_to_fmt(snv_in_fname.c_str(), NULL);
		}
	}

	// Output SNV file options
	std::string snv_out_fname = default_snv_out_fname; // "out.snv"
	snv::FMTFLAGS_T snv_out_fmt = snv::F_NULL;
	if (options[OUT_SNVFILE]) {
		snv_out_fname = options[OUT_SNVFILE].arg;
		// opts_to_fmt deduces format either from an explicitly set option or from the filename
		if (options[OUT_SNVFTYPE]) {
			snv_out_fmt = hts::snv::opts_to_fmt(snv_out_fname.c_str(), options[OUT_SNVFTYPE].arg);
		} else {
			snv_out_fmt = hts::snv::opts_to_fmt(snv_out_fname.c_str(), NULL);
		}
	} else if (options[OUT_SNVFTYPE]) {
		snv_out_fmt = hts::snv::opts_to_fmt(NULL, options[OUT_SNVFTYPE].arg);
		snv_out_fname += fmt_to_xtn(snv_out_fmt);
	}

	// Statistical model options, default is Bayesian
	ModelType model = default_model; // FREQ
	if(options[MODEL]) {
		if (strcmpi(options[MODEL].arg, "bayes") == 0) {
			model = BAYES;
		} else if (strcmpi(options[MODEL].arg, "freq") == 0) {
			model = FREQ;
		} else {
			cerr << "Error: " << options[MODEL].arg << " is not a recognized model.\n";
			success = false;
		}
	}
	
	// For frequentist identification, whether or not to allow phi to vary during analysis
	bool fixed_phi = (bool) options[FIXED_PHI];

	// Significance level option, setting the threshold for significance for identification
	// NB: pvals must be /below/ the threshold to be significant, while posterior probs must be /above/
	double sig_level;
	if (options[SIGLEVEL]) {
		sig_level = strtod(options[SIGLEVEL].arg, NULL);
	} else {// set defaults
		if (model == BAYES) {
			sig_level = default_siglevel_bayes; // .95
		} else if (model == FREQ) {
			sig_level = default_siglevel_freq; // .05
		}
	}
	
	//
	// Command-specific argument checks
	//
	
	// Now that all options are set, check for missing required args, invalid combinations, etc.
	if (quantifying) {

		if (ref_fname.empty()) {
			std::cerr <<
				"Error: reference sequence argument (-f, --reference-file) is mandatory for damage quantification.\n";
			success = false;
		}

		// Warn about ignored arguments
		// Array and counter defined in options.hpp
		warn_ignored_args(options, ident_only_args, ident_arg_count, "damage identification");
		warn_ignored_args(options, freq_ident_only_args, freq_ident_arg_count, "damage identification");
		warn_ignored_args(options, bayes_ident_only_args, bayes_ident_arg_count, "damage identification");
		if (model == FREQ) {
			warn_ignored_args(options, bayes_quant_only_args, bayes_quant_arg_count, "the Bayesian model");
		}

	} else if (identifying) {

		if (snv_in_fname.empty()) {
			std::cerr <<
				"Error: SNV file argument (-V, --snv-file) is mandatory for damage identification.\n";
			success = false;
		} else if (snv_in_fmt == snv::F_NULL) {
			std::cerr << 
				"Error: Input SNV file type not given and could not deduce type from given filename. Use -I to specify type.\n";
			success = false;
		} else if (snv_out_fmt == snv::F_NULL) {
			std::cerr <<
				"Error: Output SNV file type not given and could not deduce type from given filename (if any). Use -O to specify type.\n";
			success = false;
		}

		// Warn about ignored arguments
		// Array and counter defined in options.hpp, function defined in print_help.hpp
		warn_ignored_args(options, quant_only_args, quant_arg_count, "damage quantification");
		warn_ignored_args(options, bayes_quant_only_args, bayes_quant_arg_count, "damage quantification");
		if (model == FREQ) {
			warn_ignored_args(options, bayes_ident_only_args, bayes_ident_arg_count, "the Bayesian model");
		} else if (model == BAYES) {
			warn_ignored_args(options, freq_ident_only_args, freq_ident_arg_count, "the frequentist model");
		}
	}

	// Numeric parameter flags
	double phi = default_phi; // .01
	double alpha = default_alpha; // 1.0
	double beta = default_beta; // 1.0
	double prior_alt = default_prior_alt; // 0.5
	if (model == BAYES) {
		if (identifying) {
			// alpha hyperparameter
			if (options[ALPHA]) {
				alpha = strtod(options[ALPHA].arg, NULL);
			} else {
				global_log.v(1) <<
					"Warning: no estimate of alpha_phi was supplied. Default value of " << alpha << " will be used.\n";
			}
			// beta hyperparameter
			if (options[BETA]) {
				beta = strtod(options[BETA].arg, NULL);
			} else {
				global_log.v(1) <<
					"Warning: no estimate of beta_phi was supplied. Default value of " << beta << " will be used.\n";
			}
			// alt prior prob
			if (options[ALTPRI]) {
				prior_alt = strtod(options[ALTPRI].arg, NULL);
			} else {
				global_log.v(1) <<
					"Warning: no alternative prior probability was supplied. Default value of " << prior_alt << " will be used.\n";
			}
		}
	} else if (model == FREQ) {
		if (identifying) {
			// phi estimate
			if (options[PHI]) {
				phi = strtod(options[PHI].arg, NULL);
			} else {
				global_log.v(1) <<
					"Warning: no estimate of phi was supplied. Default value of " << phi << " will be used.\n";
			}
		}
	}
	
	// Misc options that should be sorted
	// Minimum mapping quality for analyzed reads
	int min_mapq = default_min_mapq; // 30
	if (options[MIN_MAPQ]) {
		min_mapq = atoi(options[MIN_MAPQ].arg);
	}
	// Minimum base quality for analyzed reads
	int min_baseq = default_min_baseq; // 20
	if (options[MIN_BASEQ]) {
		min_baseq = atoi(options[MIN_BASEQ].arg);
	}
	// Maximum number of reads to analyze during quantification
	// For Bayesian quant, may need to be increased to achieve convergence
	long int max_reads = default_max_reads; // 1000000
	if (options[MAX_READS]) {
		max_reads = atol(options[MAX_READS].arg);
	}
	// Symmetric log-space bounds for minimizing in freq ident process
	int minz_bound = default_minz_bound; // 15
	if (options[MINZ_BOUND]) {
		minz_bound = atoi(options[MINZ_BOUND].arg);
	}
	// Epsilon for convergence
	double eps = default_eps; // 1e-6
	if (options[EPS]) {
		eps = strtod(options[EPS].arg, NULL);
	}

	// misc flags
	// Keep duplicate reads?
	bool keep_dup = default_keep_dup; // false
	if (options[KEEP_DUP]) {
		keep_dup = options[KEEP_DUP].last()->type() == t_ON;
	}

	// error on unknown options
	if (options[UNKNOWN]) {
		success = false;
		Option* opt = options[UNKNOWN].first();
		while (opt) {
			std::cerr << "Error: Option " << opt->name << " is not recognized.\n";
			opt = opt->next();
		}
	}

	// Exit if a fatal error was encountered, after providing some help
	if (!success) {
		print_help(NULL);
		return 1;
	}
	
	//
	// Direct the program according to parsed options
	//
	using namespace hts;
	
	//
	// Damage quantification block
	//
	if (quantifying) {
		// open BAM data file
		piler p;
		if (!p.open(align_fname.c_str())) {
			std::cerr <<
				"Error: could not open BAM file \'" << align_fname << "\'.\n";
			return 1;
		}

		// set query filter parameters
		// minimum mapping and base quality for inclusion
		p.qfilter.min_mapq = min_mapq;
		p.qfilter.min_baseq = min_baseq;
		// include only reads which are properly aligned
		p.qfilter.enable_prereq_flags(BAM_FPROPER_PAIR); 
		// exclude any reads whose mates are unmapped or are supplementary reads
		p.qfilter.enable_excl_flags(BAM_FMUNMAP | BAM_FSUPPLEMENTARY); 
		// exclude duplicate reads, depending on keep_dup parameter
		if (keep_dup) {
			p.qfilter.disable_excl_flags(BAM_FDUP);
		}
		// exclude any reads which map to the same strand as their mate
		p.qfilter.excl_tandem_reads = true;
		// minimum and maximum insert size for inclusion
		p.qfilter.check_isize = true;
		p.qfilter.min_isize = 60;
		p.qfilter.max_isize = 600;

		// allocate some space for the read buffer in the piler (not a cap)
		p.reserve(100);

		// open reference sequence file
		faidx_reader faidx;
		if (!faidx.open(ref_fname.c_str())) {
			std::cerr <<
				"Error: could not open reference sequence file \'" << ref_fname << "\'.\n";
			return 1;
		}

		// initialize JSON output object
		simpleson::json::jobject quant_results;
		quant_results["bam_file"] = align_fname;
		quant_results["damage_type"] = dtype;

		if (model == FREQ) {
			if (!plain_output) {
				std::cerr << "Starting frequentist quantification...\n";
			}
			// do quantification (-->frontend/orient-bias-quantify.hpp)
			success = orient_bias_quantify_freq(ref, alt, p, faidx, quant_results, max_reads, plain_output);

			if (!success) {
				std::cerr << "Quantification process failed.\n";
				return 1;
			}
		} else if (model == BAYES) {
			if (!plain_output) {
				std::cerr << "Starting Bayesian quantification (this will take a bit)...\n";
			}
			// do quantification (-->frontend/orient-bias-quantify.hpp)
			if (options[ALPHA] && options[BETA]) {
				// use initial estimates for alpha and beta if they were passed (both reqd)
				success = orient_bias_quantify_bayes(ref, alt, p, faidx, quant_results, max_reads, plain_output, eps, alpha, beta);
			} else {
				if (options[ALPHA] || options[BETA]) {
					frontend::global_log.v(1) << "Warning: Both --alpha and --beta must be specified to pass an initial estimate to Bayesian quant.\n";
				}
				// if no estimates given, use frequentist quant to obtain an initial estimate
				success = orient_bias_quantify_bayes(ref, alt, p, faidx, quant_results, max_reads, plain_output, eps);
			}
				
			if (!success) {
				std::cerr << "Quantification process failed.\n";
				return 1;
			}
		}

		// write out to JSON if asked to
		if (!json_fname.empty()) {
			string json_out = (string) quant_results;
			indent_serialized_json(json_out); // defined in frontend/string.hpp
			ofstream writeout (json_fname.c_str());
			writeout << json_out << endl;
			writeout.close();
		}
	} // End quantification block

	//
	// Damage identification block
	// 
	if (identifying) {
		if (!plain_output) {
			global_log.v(1) << "Starting identification...\n";
		}
		
		// open BAM data file
		fetcher alignment_file;
		if (!alignment_file.open(align_fname.c_str())) {
			std::cerr <<
				"Error: Could not open BAM file \'" << align_fname << "\'.";
			return 1;
		}
		
		// streamer class helps pass paramters from input SNV to output SNV file
		snv::streamer snv_files (snv_in_fname.c_str(), snv_out_fname.c_str(), snv_in_fmt, snv_out_fmt);
		
		// do the identification (-->frontend/orient-bias-identify.hpp)
		if (model == BAYES) {
			success = orient_bias_identify_bayes(ref, alt, eps, minz_bound, alpha, beta,
				prior_alt, sig_level, alignment_file, snv_files, plain_output);
		} else if (model == FREQ) {
			success = orient_bias_identify_freq(ref, alt, eps, minz_bound, phi,
				sig_level, alignment_file, snv_files, fixed_phi, plain_output);
		}
		
		if (!success) {
			global_log.v(1) << "Identification process failed.";
			return 1;
		}
	} // identification block

	return 0;

}// main
