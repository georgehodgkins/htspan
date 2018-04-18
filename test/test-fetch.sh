#!/bin/bash

bam=../data/test.bam
target=chr17
bindir=../bin
outdir=out/fetch
testdir=ans/fetch

poss=( 7674420 7674360 7674361 )

mkdir -p $outdir

for pos in ${poss[@]}; do
	$bindir/hts-fetch $bam $target $pos 0 0 0 > ${outdir}/out_${target}-${pos}-nodup.txt
	$bindir//hts-fetch $bam $target $pos 1 0 0 > ${outdir}/out_${target}-${pos}-dup.txt
done


echo "Testing diff ..."

dups=( dup nodup )

for dup in ${dups[@]}; do
	for pos in ${poss[@]}; do
		diff ${outdir}/out_${target}-${pos}-${dup}.txt \
			${testdir}/test_${target}-${pos}-${dup}.txt
		echo "${outdir}/out_${target}-${pos}-${dup}.txt ... $?"
	done
done
