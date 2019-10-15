#!/usr/bin/env Rscript

library(io)
library(argparser)

pr <- arg_parser("Convert GATK orientation bias VCF to SNV")
pr <- add_argument(pr, "input", help="input VCF file");
pr <- add_argument(pr, "--output", help="output tab-delimited file");

argv <- parse_args(pr);

if (is.na(argv$output)) {
	output <- set_fext(as.filename(argv$input), ext="snv");
} else {
	output <- argv$output;
}

input <- argv$input;

x <- qread(input);

d <- x[, c(1, 2, 4, 5)];
d$OBP <- unlist(lapply(x$ffpe,
	function(z) {
		if (is.null(z$OBP)) {
				NA
		} else {
			z$OBP[1]
		}
	}
));

qwrite(d, output, type="tsv");

