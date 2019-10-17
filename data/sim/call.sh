#!/bin/bash

# Call variants on aligned reads.

set -euo pipefail

if (( $# < 2 )); then
	echo "usage : $0 <bam> <reference-fasta> [mut_rate]" >&2
	exit 1
fi

home="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

bam=$1
ref=$2

max_depth=5000
mut_rate=${3:-0.99}

prefix=${bam%.*}

samtools mpileup -Ou -I -d $max_depth -f $ref $bam |
	bcftools call -mv -Ov -P $mut_rate -o $prefix.calls.vcf
$home/vcf2snv.sh $prefix.calls.vcf $prefix.calls.snv

# consensus-caller
#	bcftools call -cv -Ov -p 1 -o $prefix.calls.vcf

