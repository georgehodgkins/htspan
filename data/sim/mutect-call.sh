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

prefix=${bam%.*}

if [[ ! -f ${ref%.*}.dict ]]; then
	gatk CreateSequenceDictionary -R $ref -O ${ref%.*}.dict
fi

gatk Mutect2 -R $ref -I $bam -O $prefix.raw.calls.vcf
gatk FilterMutectCalls -R $ref -V $prefix.raw.calls.vcf -O $prefix.ann.calls.vcf
gatk SelectVariants --exclude-filtered -V $prefix.ann.calls.vcf -O $prefix.calls.vcf
rm $prefix.raw.calls.vcf $prefix.ann.calls.vcf

$home/vcf2snv.sh $prefix.calls.vcf $prefix.calls.snv

