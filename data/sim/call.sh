#!/bin/bash

# Call variants on aligned reads.

set -euo pipefail

if (( $# < 2 )); then
	echo "usage : $0 <bam> <reference-fasta>" >&2
	exit 1
fi

bam=$1
ref=$2

mut_rate=${3:-0.02}

fname=${bam##*/}
prefix=${fname%.*}

bcftools mpileup -Ou -f $ref $bam |
	bcftools call -mv -Ov -P $mut_rate -o $prefix.calls.vcf
./vcf2snv.sh $prefix.calls.vcf $prefix.calls.snv

