#!/bin/bash

bam=../data/test.bam
bindir=../bin
outdir=out/pileup
testdir=ans/pileup

mkdir -p $outdir

# do not keep duplicate
$bindir/hts-pileup $bam 0 30 20 10 > ${outdir}/pileup.tsv


echo "Testing diff (hts-pileup) ..."

diff ${outdir}/pileup.tsv ${testdir}/pileup.tsv
echo "${outdir}/pileup.tsv ... $?"

