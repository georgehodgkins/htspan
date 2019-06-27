#!/bin/bash

# Simulate reads.

set -euo pipefail

if (( $# < 2 )); then
	echo "usage : $0 <out-prefix> <reference-fasta>" >&2
	exit 1
fi

out=$1
ref=$2

seed=${3:-0}

r1=$out.r1.fq
r2=$out.r2.fq

# 2*100 bp reads * 10000 reads / 19149 bp seq = 104x coverage
# SNV mutation rate of 0.05
# no indels
printf 'chrom\tpos\tref\talt\n' > $out.snv
wgsim -h -1 100 -2 100 -S $seed -e 0.05 -N 10000 -R 0 \
	$ref $r1 $r2 | 
	cut -f 1-4 >> $out.snv

