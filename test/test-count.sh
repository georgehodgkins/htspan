#!/bin/bash

bam=../data/test.bam
bindir=../bin
outdir=out/count
testdir=ans/count

mkdir -p $outdir

# do not keep duplicate
$bindir/hts-count $testdir/snvs.tsv $bam $outdir/cov.tsv 0 0 0 1


echo "Testing diff (hts-count) ..."

diff ${outdir}/cov.tsv ${testdir}/cov.tsv
echo "${outdir}/cov.tsv ... $?"

