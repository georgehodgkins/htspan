#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

#include "frontend/options.hpp"
#include "frontend/print-help.hpp"
#include "frontend/orient-bias-identify.hpp"
#include "frontend/orient-bias-quantify.hpp"
#include "frontend/optionparser.hpp"
#include "frontend/cstring.hpp"
#include "frontend/simul_writer.hpp"

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

// TODO: expose simulation capabilities

int main (int argc, char** argv) {
	// enum for model options
	enum ModelType {BAYES, FREQ};

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
			strcmpi(command, "-?") == 0) {
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
		print_help(options[HELP].arg);
		return 0;
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
	// TODO: update so '-' redirects to stdout
	bool use_stdout = true;
	if (options[STDOUT]) {
		use_stdout = options[STDOUT].last()->type() == t_ON;
	}
	int verbosity = 1;
	if (options[VERBOSITY]) {
		verbosity = atoi(options[VERBOSITY].arg);
	}
	bool log_to_file = (bool) options[LOGFILE];
	std::string log_fname;
	if (log_to_file) {
		log_fname = options[LOGFILE].arg;
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
	if (options[DAMAGE_TYPE]) {
		if (strcmpi(options[DAMAGE_TYPE].arg, "ffpe") == 0) {
			ref = nuc_C;
			alt = nuc_T;
		} else if (strcmpi(options[DAMAGE_TYPE].arg, "oxog") == 0) {
			ref = nuc_G;
			alt = nuc_T;
		}
	} else if (options[REF] && options[ALT]) {
		ref = char_to_nuc(options[REF].arg[0]);
		alt = char_to_nuc(options[ALT].arg[0]);
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

	// Simulation options
	bool internal_sim = (bool) options[INT_SIM];
	bool external_sim = (bool) options[EXT_SIM];
	std::string ext_sim_fname;
	if (external_sim) {
		ext_sim_fname = options[EXT_SIM].arg;
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
	// TODO: add extension if a name is given as well?
	std::string snv_out_fname = "out"; // default name
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
	ModelType model = BAYES;
	if(options[MODEL]) {
		if (strcmpi(options[MODEL].arg, "freq") == 0) {
			model = FREQ;
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
			sig_level = .95;
		} else if (model == FREQ) {
			sig_level = .05;
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
		for (size_t n = 0; n < ident_arg_count; ++n) {
			if (options[ident_only_args[n]]) {
				global_log.v(1) << "Warning: option [-/--]" << options[ident_only_args[n]].name << " only applies to damage identification. Ignored.\n";
			}
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
		// Array and counter defined in options.hpp
		for (size_t n = 0; n < quant_arg_count; ++n) {
			if (options[quant_only_args[n]]) {
				global_log.v(1) << "Warning: option [-/--]" << options[quant_only_args[n]].name << " only applies to damage quantification. Ignored.\n";
			}
		}
	}

	// Numeric parameter flags
	// TODO: use --alpha/beta for initial estimates for quant as well
	double phi = .01;
	double alpha = 1.0;
	double beta = 1.0;
	double prior_alt = .5;
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
	int min_mapq = 30;
	if (options[MIN_MAPQ]) {
		min_mapq = atoi(options[MIN_MAPQ].arg);
	}
	// Minimum base quality for analyzed reads
	int min_baseq = 20;
	if (options[MIN_BASEQ]) {
		min_baseq = atoi(options[MIN_BASEQ].arg);
	}
	// Maximum number of reads to analyze during quantification
	// For Bayesian quant, may need to be increased to achieve convergence
	long int max_reads = 1000000;
	if (options[MAX_READS]) {
		max_reads = atol(options[MAX_READS].arg);
	}
	// Symmetric log-space bounds for minimizing in freq ident process
	int minz_bound = 15;
	if (options[MINZ_BOUND]) {
		minz_bound = atoi(options[MINZ_BOUND].arg);
	}
	// Epsilon for convergence
	// TODO: apply this to Bayesian quant?
	double eps = 1e-6;
	if (options[EPS]) {
		eps = strtod(options[EPS].arg, NULL);
	}

	// simulation-only flags
	double theta_sim = 0.1;
	double phi_sim = 0.1;
	double err_mean_sim = 30.0;
	double err_sd_sim = 2.0;
	if (internal_sim) {
		if (external_sim) {
			global_log.v(1) << "Warning: the -S/--internal-sim flag overrides -s/--external-sim.\n";
		}
		if (!options[THETA_SIM] || !options[PHI_SIM]
			|| !options[ERR_MEAN_SIM] || !options[ERR_SD_SIM]) {
			std::cerr <<
				"Error: All --*-sim parameters are required for internal simulation. Use --help for usage.\n";
			success = false;
		} else {
			theta_sim = strtod(options[THETA_SIM].arg, NULL);
			phi_sim = strtod(options[PHI_SIM].arg, NULL);
			err_mean_sim = strtod(options[ERR_MEAN_SIM].arg, NULL);
			err_sd_sim = strtod(options[ERR_SD_SIM].arg, NULL);
		}
	} else {
		for (size_t n = 0; n < sim_arg_count; ++n) {
			if (options[intsim_only_args[n]]) {
				global_log.v(1) << "Warning: option [-/--]" << options[intsim_only_args[n]].name << " only applies to internal simulation. Ignored.\n";
			}
		}
	}

	// misc flags
	// Keep duplicate reads?
	bool keep_dup = false;
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
		if (internal_sim) {
			global_log.v(1) << "Warning: Quantification internal sim code does not exist yet.\n";
			return 0;
		} else if (external_sim) {
			global_log.v(1) << "Warning: Quantification external sim code does not exist yet.\n"; 
			return 0;
		} else {
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
			p.qfilter.min_isize = 60;
			p.qfilter.max_isize = 600;

			// open reference sequence file
			faidx_reader faidx;
			if (!faidx.open(ref_fname.c_str())) {
				std::cerr <<
					"Error: could not open reference sequence file \'" << ref_fname << "\'.\n";
				return 1;
			}

			if (model == FREQ) {
				if (!plain_output) {
					std::cerr << "Starting frequentist quantification...\n";
				}
				// do quantification (-->frontend/orient-bias-quantify.hpp)
				success = orient_bias_quantify_freq(ref, alt, min_mapq, min_baseq,
					keep_dup, max_qreads, p, faidx, plain_output);

				if (!success) {
					std::cerr << "Quantification process failed.\n";
					return 1;
				}
			} else if (model == BAYES) {
				if (!plain_output) {
					std::cerr << "Starting Bayesian quantification (this will take a bit)...\n";
				}
				// do quantification (-->frontend/orient-bias-quantify.hpp)
				success = orient_bias_quantify_bayes(ref, alt, min_mapq, min_baseq,
					keep_dup, max_qreads, p, faidx, plain_output);
					
				if (!success) {
					std::cerr << "Quantification process failed.\n";
					return 1;
				}
			}
		}
	} // End quantification block
	//
	// Damage identification block
	// 
	if (identifying) {
		if (internal_sim) { // TODO: move simulation calls to intermediate file
			//orient_bias_filter_f obfilter(alt, ref, 0);
			//obfilter.simulate(theta_sim, phi_sim, err_mean_sim, err_sd_sim);
		} else if (external_sim) {
			//orient_bias_filter_f obfilter(alt, ref, 0,
			//	-1*minz_bound, minz_bound, minz_eps, minz_iter);// change count param when param options are added
			//obfilter.read(ext_sim_fname.c_str());
		} else {
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
		}
	} // identification block
	return 0;
}// main
