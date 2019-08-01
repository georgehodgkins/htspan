#!/bin/bash

# Call variants on aligned reads.

set -euo pipefail

if (( $# < 2 )); then
	echo "usage : $0 <bam> <reference-fasta>" >&2
	exit 1
fi

home="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

bam=$1
ref=$2

min_vaf=0.01

prefix=${bam%.*}

# use params to detect rare variants
# ensure that mnvs are split into individual snvs
freebayes -O -K -F $min_vaf -f $ref $bam |
	vcfallelicprimitives > $prefix.calls.vcf
$home/vcf2snv.sh $prefix.calls.vcf $prefix.calls.snv

