#!/bin/bash

# Simulate reads.

set -euo pipefail

if (( $# < 2 )); then
	echo "usage : $0 <out-prefix> <reference-fasta> [<seed> <mut_rate> <base_error> <rd_count>]" >&2
	exit 1
fi

out=$1
ref=$2

seed=${3:--1}
mut_rate=${4:-0.02}
base_error=${5:-0.005}
rd_count=${6:-10000}


# 2*100 bp reads * 10000 reads / 19149 bp seq = 104x coverage
# SNV mutation rate of 0.05
# no indels
printf 'chrom\tpos\tref\talt\n' > $out.snv
dwgsim -H -1 100 -2 100 -z $seed -r ${mut_rate} -e ${base_error} -N ${rd_count} -R 0 \
	$ref $out
cat $out.mutations.txt | cut -f 1-4 >> $out.snv

# the rest of the scripts expect these names
mv $out.bwa.read1.fastq $out.r1.fq
mv $out.bwa.read2.fastq $out.r2.fq

# don't need these files
rm -f $out.mutations.txt $out.mutations.vcf
