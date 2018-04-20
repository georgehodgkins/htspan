#!/bin/bash

bam=../data/test.bam
snv=../data/snv.tsv
target=chr17
bindir=../bin
outdir=out/orient-bias-filter
testdir=ans/orient-bias-filter

mkdir -p $outdir

$bindir/hts-orient-bias-filter $bam $snv 0.01 0 0 0 > ${outdir}/out.vtr


echo "TODO Testing diff ... "
cat ${outdir}/out.vtr
