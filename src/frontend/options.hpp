#ifndef _HTSPAN_OPTION_HPP_
#define _HTSPAN_OPTION_HPP_

#include "optionparser.hpp"
#include "arg.hpp"

/**
* This file contains the list of command-line flags
* that can be passed to hts-orient-bias.
*
* It should be mostly self-documenting, since each option
* includes a help text describing its usage.
*/

namespace hts {
	
namespace frontend {

// These indices uniquely identify options
enum OptionIndex {UNKNOWN=0, REF=1, ALT=2,VERBOSITY=3, LOGFILE=4, BAMFILE=5,
	REFFILE=6, IN_SNVFILE=7, OUT_SNVFILE=8, PHI=9, STDOUT=10,
	MIN_MAPQ=11, MIN_BASEQ=12, KEEP_DUP=13, MAX_READS=14, MINZ_BOUND=15, EPS=16,
	DAMAGE_TYPE=17, ALPHA=18, BETA=19, ALTPRI=20, MODEL=21, IN_SNVFTYPE=22, OUT_SNVFTYPE=23,
	SIGLEVEL=24, FIXED_PHI=25, PLAIN=26, JSON_OUT=27, HELP=28};

enum OptionType {t_ON, t_OFF, t_OTHER};

// These arrays track arguments that apply to both commands, for printing help texts
// common_args help will be printed before args unique to the command
const OptionIndex common_args[] = {DAMAGE_TYPE, BAMFILE, REF, ALT, MODEL};
const size_t common_arg_count = sizeof(common_args)/sizeof(OptionIndex);

const OptionIndex utility_args[] = {VERBOSITY, LOGFILE, STDOUT, PLAIN, HELP};
const size_t utility_arg_count = sizeof(utility_args)/sizeof(OptionIndex);

const OptionIndex filter_args[] = {KEEP_DUP, MIN_MAPQ, MIN_BASEQ};
const size_t filter_arg_count = sizeof(filter_args)/sizeof(OptionIndex);

// These arrays track arguments that only apply to bayes/freq quantification/identification
// Used to report ignored arguments to the user and for printing help texts
const OptionIndex quant_only_args[] = {REFFILE, MAX_READS, JSON_OUT};
const size_t quant_arg_count = sizeof(quant_only_args)/sizeof(OptionIndex);

const OptionIndex bayes_quant_only_args[] = {ALPHA, BETA, EPS};
const size_t bayes_quant_arg_count = sizeof(bayes_quant_only_args)/sizeof(OptionIndex);

const OptionIndex ident_only_args[] = {IN_SNVFILE, OUT_SNVFILE, IN_SNVFTYPE, OUT_SNVFTYPE, MINZ_BOUND, SIGLEVEL, EPS};
const size_t ident_arg_count = sizeof(ident_only_args)/sizeof(OptionIndex);

const OptionIndex freq_ident_only_args[] = {PHI, FIXED_PHI};
const size_t freq_ident_arg_count = sizeof(freq_ident_only_args)/sizeof(OptionIndex);

const OptionIndex bayes_ident_only_args[] = {ALPHA, BETA, ALTPRI};
const size_t bayes_ident_arg_count = sizeof(bayes_ident_only_args)/sizeof(OptionIndex);

/**
* Main array of options passed to the parser.
* This array must stay ordered least to greatest index for help printing to work.
* Note that in help descrs, \r feeds a newline + indentation for all subsequent lines in that description
* All option checks are defined in arg.hpp, except for built-in Arg::None and Arg::Optional.
*/
const option::Descriptor usage[] = {
	{
		UNKNOWN,
		t_OTHER,
		"",
		"",
		Arg::None,
		"htspan: Corrects errors in high-throughput sequencing data."
	},{
		REF,
		t_OTHER, 
		"R",
		"reference-nuc",
		Arg::CanonicalNucleotide,
		"-R, --reference-nuc\rThe reference nucleotide on read 1 (A, C, G, or T).\n"
		"This and -A will be overridden by the -t flag."
	},{
		ALT,
		t_OTHER,
		"A",
		"alternate-nuc",
		Arg::CanonicalNucleotide,
		"-A, -alternate-nuc\rThe alternate nucleotide on read 1 (A, C, G, or T).\n"
		"This and -R will be overridden by the -t flag."
	},{
		VERBOSITY,
		t_OTHER,
		"v",
		"verbosity",
		Arg::Verbosity,
		"-v, --verbosity [1]\rLevel of runtime verbosity (separate from result output), "
		"from 0 (silent) - 3 (details about every pushed read).\n"
		"NB: Level 3 may perform extra computations to produce its output, "
		"at a penalty to overall program performance."
	},{
		LOGFILE,
		t_OTHER,
		"l",
		"log-file",
		Arg::OutputFile,
		"-l, --log-file\rPath to write runtime output to.\n"
		"Note: Runtime output and results on stdout are controlled by the -O flag."
	},{
		BAMFILE,
		t_OTHER,
		"b",
		"alignment-file",
		Arg::PairedEndAlignmentFile,
		"-b, --alignment-file\rPath for the BAM data being examined for damage. "
		"Must be paired-end."
	},{
		REFFILE,
		t_OTHER,
		"f",
		"reference-file",
		Arg::ReferenceFile,
		"-f, --reference-file\rPath for the reference genome (FASTA or FASTQ) "
		"to compare against in damage quantification."
	},{
		IN_SNVFILE,
		t_OTHER,
		"V",
		"in-snv-file",
		Arg::InSnvFile,
		"-V, --snv-file\rPath for the list of SNVs to be examined in damage identification.\n"
		"Valid formats are plain TSV or VCF/BCF (can be compressed with gzip or bgzip).\n"
		"Valid extensions are:\'snv\', \'tsv\', \'vcf\', \'bcf\', \'gz\', or \'bgz\'."
	},{
		OUT_SNVFILE,
		t_OTHER,
		"o",
		"out-snv-file",
		Arg::OutputFile,
		"-o, --out-snv-file [out.snv]\rPath for the output list of filtered SNVs.\n"
		"If the -O flag is not set, output file type will be deduced from the extension of this file.\n"
	},{
		PHI,
		t_OTHER,
		"p",
		"phi",
		Arg::Probability,
		"-p, --phi [.01]\rEstimate of global damage for frequentist damage identification."
	},{
		STDOUT,
		t_ON,
		"",
		"stdout",
		Arg::None,
		"--stdout [true]\rTurns on printing runtime output and results to stdout, "
		"in addition to any log or result files specified "
		"(default behavior)."
	},{
		STDOUT,
		t_OFF,
		"",
		"no-stdout",
		Arg::None,
		"--no-stdout [false]\rTurns off printing runtime output and results to stdout.\n"
		"Note that if you use this flag and do not specify a result file (-o), "
		"you will not be able to see your results."
	},{
		MIN_MAPQ,
		t_OTHER,
		"",
		"min-mapq",
		Arg::MinMapQ, 
		"--min-mapq [30]\rMinimum mapping quality for a read to be included in analysis."
	},{
		MIN_BASEQ,
		t_OTHER,
		"",
		"min-baseq",
		Arg::MinBaseQ, 
		"--min-baseq [20]\rMinimum base quality for a read to be included in analysis."
	},{
		KEEP_DUP,
		t_ON,
		"D",
		"keep-dup",
		Arg::None,
		"-D, --keep-dup [false]\rKeep reads at the same reference position present in multiple sequences."
	},{
		MAX_READS,
		t_OTHER,
		"",
		"max-reads",
		Arg::MaxQReads, 
		"--max-reads [1e6]\rMaximum number of reads to process in the damage quantification process. "
		"Note that this count only includes reads that pass the quality filter and are relevant to the analysis."
	},{
		MINZ_BOUND,
		t_OTHER,
		"",
		"minimizer-bound",
		Arg::MinzBound, 
		"--minimizer-bound [15]\rMagnitude (log space) of the symmetrical domain to minimize in during damage identification. "
		"i.e. for the default arg the minimizer will search from xval e^-15 to e^15 in constrained space."
	},{
		EPS,
		t_OTHER,
		"e",
		"epsilon",
		Arg::Eps,
		"-e, --epsilon [1e-6]\rThreshold for convergence in optimization and integration methods used."
	},{
		DAMAGE_TYPE,
		t_OTHER,
		"t",
		"damage-type",
		Arg::DamageType,
		"-t, --damage-type\rType of damage to identify or quantify (e.g. what variant type to analyze). "
		"Choices are \'ffpe\' or \'oxog\'. This flag will override variants selected using -R and -A."
	},{
		ALPHA,
		t_OTHER,
		"",
		"alpha",
		Arg::PositiveDouble,
		"--alpha [1]\rAlpha hyperparameter for beta distribution; used as a parameter in Bayesian identification "
		"and as an initial estimate in Bayesian quantification.\nNote that both --alpha and --beta must be set to be used as initial estimates."
	},{
		BETA,
		t_OTHER,
		"",
		"beta",
		Arg::PositiveDouble,
		"--beta [1]\rBeta hyperparameter for beta distribution; used as a parameter in Bayesian identification "
		"and as an initial estimate in Bayesian quantification.\nNote that both --alpha and --beta must be set to be used as initial estimates."
	},{
		ALTPRI,
		t_OTHER,
		"",
		"alt-prior",
		Arg::Probability,
		"--alt-prior [.5]\rPrior probability of the alternative hypothesis (theta != 0) in the Bayesian identification model."
	},{
		MODEL,
		t_OTHER,
		"M",
		"model",
		Arg::Model,
		"-M, --model [freq]\rModel to use for identification or quantification. "
		"Choices are \'freq\' or \'bayes\'."
	},{
		IN_SNVFTYPE,
		t_OTHER,
		"I",
		"in-snv-type",
		Arg::SnvFType,
		"-I, --in-snv-type [auto]\rFormat of SNV input file selecting variants for damage identification.\n"
		"Choices are \'tsv\' or \'vcf\'. If this option is not provided, type will be deduced from the file extension.\n"
		"For input files, compression is supported and autodetected.\n"
		"NB: The \'vcf\' argument also supports BCF files."
	},{
		OUT_SNVFTYPE,
		t_OTHER,
		"O",
		"out-snv-type",
		Arg::SnvFType,
		"-O, --out-snv-type [auto]\rFormat of SNV output file containing filtered variants.\n"
		"Choices are \'tsv\', \'vcf\', or \'vcf-bgz\'. If this option is not provided, type will be deduced from the file extension, "
		"with the .gz and .bgz extensions both indicating BGZF compression."
	},{
		SIGLEVEL,
		t_OTHER,
		"g",
		"sig-level",
		Arg::Probability,
		"-g, --sig-level [.05 freq/.95 bayes]\rThreshold for a variant to be considered damaged by the filter.\n"
		"Specifically, a variant with a p-val /above/ (freq model) or a posterior probability /below/ (bayes model) "
		"this threshold is flagged as damaged by the filter.\n"
		"Set to 0 to disable automatic filtering (statistics will still be recorded in the output file)."
	},{
		FIXED_PHI,
		t_OTHER,
		"",
		"fixed-phi",
		Arg::None,
		"--fixed-phi\rIn the frequentist model, fix phi to the supplied value "
		"(or the default if none was supplied) rather than estimating it."
	},{
		PLAIN,
		t_OTHER,
		"",
		"plain",
		Arg::None,
		"--plain\rMake output more machine-readable"
	},{
		JSON_OUT,
		t_OTHER,
		"J",
		"json-out",
		Arg::OutputFile,
		"-J, --json-out\rOutput quantification results and extra info to a JSON file."
	},{
		HELP,
		t_OTHER,
		"?",
		"help",
		Arg::Optional,
		"-?, --help\rPrint a help message (pass a command, argument name, or \"utility\" to get more detailed help)."
	},{// null terminator
		UNKNOWN,
		0,
		0,
		0,
		0,
		0
	}
};// option array

size_t total_arg_count = sizeof(usage)/sizeof(option::Descriptor);

}// namespace frontend

}// namespace htspan

#endif // _HTSPAN_OPTION_HPP_
