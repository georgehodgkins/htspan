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

#include "htspan/nucleotide.hpp"
#include "htspan/fetcher.hpp"
#include "htspan/piler.hpp"
#include "htspan/print.hpp"

#include "htspan/io/faidx_reader.hpp"
#include "htspan/io/snv_reader.hpp"
#include "htspan/io/simul_writer.hpp"

using namespace hts::frontend;

// TODO Update classes to use proper logging and verbosity
// TODO Eliminate string/cstr juggling

int main (int argc, char** argv) {
	// used for clarity when setting model and integrator options
	enum ModelType {BAYES, FREQ};
	enum SnvFileType {TSV, VCF};

	//
	// Process input arguments
	// 
	if (argc <= 1) { // i.e. no arguments given
		print_help(NULL);
		return 1;
	}
	argv++; // ignore executable name
	argc--;
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
	// so any args here are assumed valid

	//
	// Argument parsing block
	// 
	bool success = true;

	// Set up logging and result output
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
	bool results_to_file = (bool) options[RESFILE];
	std::string result_fname;
	if (results_to_file) {
		result_fname = options[RESFILE].arg;
	}

	// Direct file and result output appropriately
	// Note that fatal errors will go to stdout regardless either by
	// cerr in the frontend or thrown exceptions in the backend
	if (log_to_file) {
		global_log.add_file(log_fname);
	}
	global_log.use_cerr(use_stdout);
	// Verbosity levels: 0=silent, 1=warnings only 2=some runtime info 3=too much runtime info
	global_log.set_verbosity(verbosity);
	if (results_to_file) {
		global_out.add_file(result_fname);
	}
	global_out.use_cout(use_stdout);
	// disable verbosity on result output
	global_out.set_verbosity(-1);
	if (!use_stdout && !results_to_file) {
		global_log <<
		"Warning: There is no set output method for the results; results will not be accessible.\n";
	}

	// Set reference and alternative nucleotides and check that they differ
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

	// Simulation flags
	bool internal_sim = (bool) options[INT_SIM];
	bool external_sim = (bool) options[EXT_SIM];
	std::string ext_sim_fname;
	if (external_sim) {
		ext_sim_fname = options[EXT_SIM].arg;
	}

	// Input file flags
	std::string align_fname;
	if (!options[BAMFILE]) {
		std::cerr << 
			"Error: alignment file argument (-b, --alignment-file) is required.\n";
		success = false;
	} else {
		align_fname = options[BAMFILE].arg;
	}
	std::string ref_fname;
	if (options[REFFILE]) {
		ref_fname = options[REFFILE].arg;
	}
	std::string snv_fname;
	SnvFileType snv_xtn;
	if (options[SNVFILE]) {
		snv_fname = options[SNVFILE].arg;
	}

	// get statistical model (default is bayes)
	// 0 = bayes, 1 = freq
	ModelType model = BAYES;
	if(options[MODEL]) {
		if (strcmpi(options[MODEL].arg, "freq") == 0) {
			model = FREQ;
		}
	}
	

	// argument checks for specific commands
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
		if (snv_fname.empty()) {
			std::cerr <<
				"Error: SNV file argument (-V, --snv-file) is mandatory for damage identification.\n";
			success = false;
		} else {
			// check if SNV file type is set manually
			if (options[SNVFTYPE]) {
				if (strcmpi(options[SNVFTYPE].arg, "tsv") == 0) {
					snv_xtn = TSV;
				} else {// VCF or BCF
					snv_xtn = VCF;
				}
			} else { // deduce type if it is not set manually
				const char* xtn = strrchr(snv_fname.c_str(), '.');
				if (xtn != NULL) {
					if (strcmpi(&xtn[1], "tsv") == 0) {
						snv_xtn = TSV;
					} else if (strcmpi(&xtn[1], "vcf") == 0 ||
							strcmpi(&xtn[1], "bcf") == 0) {
						snv_xtn = VCF;
					} else {
						std::cerr << "Could not deduce SNV file type from extension. Please specify a type with --snv-type.";
						success = false;
					}
				} else { // no extension on filename
					std::cerr << "Could not deduce SNV file type from extension. Please specify a type with --snv-type.";
					success = false;
				}
			}
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
	double phi = .01;
	double alpha = .1;
	double beta = .1;
	double prior_alt = .5;
	if (model == BAYES) {
		// alpha hyperparameter
		if (options[ALPHA]) {
			alpha = strtod(options[ALPHA].arg, NULL);
		} else {
			global_log.v(1) <<
				"Warning: no estimate of alpha was supplied. Default value of .1 will be used.\n";
		}
		// beta hyperparameter
		if (options[BETA]) {
			beta = strtod(options[BETA].arg, NULL);
		} else {
			global_log.v(1) <<
				"Warning: no estimate of beta was supplied. Default value of .1 will be used.\n";
		}
		// alt prior prob
		if (options[ALTPRI]) {
			prior_alt = strtod(options[ALTPRI].arg, NULL);
		} else {
			global_log.v(1) <<
				"Warning: no estimate of beta was supplied. Default value of .1 will be used.\n";
		}
	} else if (model == FREQ) {
		// phi estimate
		if (options[PHI]) {
			phi = strtod(options[PHI].arg, NULL);
		} else {
			global_log.v(1) <<
				"Warning: no alternative prior probability was supplied. Default value of .5 will be used.\n";
		}
	}
	int min_mapq = 5;
	if (options[MIN_MAPQ]) {
		min_mapq = atoi(options[MIN_MAPQ].arg);
	}
	int min_baseq = 20;
	if (options[MIN_BASEQ]) {
		min_baseq = atoi(options[MIN_MAPQ].arg);
	}
	long int max_qreads = 5e7;
	if (options[MAX_QREADS]) {
		max_qreads = atol(options[MAX_QREADS].arg);
	}
	int minz_bound = 15;
	if (options[MINZ_BOUND]) {
		minz_bound = atoi(options[MINZ_BOUND].arg);
	}
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
	bool keep_dup = false;
	if (options[KEEP_DUP]) {
		keep_dup = options[KEEP_DUP].last()->type() == t_ON;
	}

	//Exit if a fatal error was encountered, after providing some help
	if (!success) {
		if (quantifying) {
			print_help("quantify");
		} else if (identifying) {
			print_help("identify");
		}
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
			global_log.v(1) << "Warning: Quantification eternal sim code does not exist yet.\n"; 
			return 0;
		} else {
			global_log.v(1) << "Starting quantification...\n";
			success = orient_bias_quantify(ref, alt, align_fname.c_str(), ref_fname.c_str(), min_mapq, min_baseq, keep_dup, max_qreads);
			if (!success) {
				global_log.v(1) << "Quantification process failed.\n";
				return 1;
			}
		}
	} // Quantification block
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
			global_log.v(1) << "Starting identification...\n";
			if (model == BAYES) {
				if (snv_xtn == TSV) {
					success = orient_bias_identify_bayes<snv::tsv_reader>(ref, alt, eps, minz_bound,
						snv_fname.c_str(), align_fname.c_str(), alpha, beta, prior_alt);
				} else if (snv_xtn == VCF) {
					success = orient_bias_identify_bayes<snv::vcf_reader>(ref, alt, eps, minz_bound,
						snv_fname.c_str(), align_fname.c_str(), alpha, beta, prior_alt);
				}
			} else if (model == FREQ) {
				if (snv_xtn == TSV) {
					success = orient_bias_identify_freq<snv::tsv_reader>(ref, alt, eps, minz_bound,
						snv_fname.c_str(), align_fname.c_str(), phi);
				} else if (snv_xtn == VCF) {
					success = orient_bias_identify_freq<snv::vcf_reader>(ref, alt, eps, minz_bound,
						snv_fname.c_str(), align_fname.c_str(), phi);
				}
			}
			if (!success) {
				global_log.v(1) << "Identification process failed.";
				return 1;
			}
		}
	} // identification block
}// main
