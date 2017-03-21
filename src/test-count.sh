#!/bin/bash

bam=../data/test.bam
outdir=out/count
testdir=test/count

mkdir -p $outdir

# do not keep duplicate
./hts-count $testdir/snvs.tsv $bam $outdir/cov.tsv 0 0 0 1


echo "Testing diff ..."

diff ${outdir}/cov.tsv ${testdir}/cov.tsv
echo "${outdir}/cov.tsv ... $?"

