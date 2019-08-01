#!/bin/bash

# Align paired-end reads

set -euo pipefail

if (( $# < 2 )); then
	echo "usage : $0 <prefix> <reference-fasta>" >&2
	exit 1
fi

prefix=$1
ref=$2

r1=$prefix.r1.fq
r2=$prefix.r2.fq


if [[ ! -f ${ref}.bwt ]]; then
	bwa index ${ref}
fi

sm=${prefix##*/}
rg="@RG\\tID:$sm\\tPU:$sm\\tSM:$sm\\tLB:$sm\\tPL:illumina"

bwa mem $ref $r1 $r2 -R $rg |
	samtools view -b > $prefix.raw.bam
samtools sort $prefix.raw.bam -o $prefix.bam
rm $prefix.raw.bam
samtools index $prefix.bam

