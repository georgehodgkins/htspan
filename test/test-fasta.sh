#!/bin/bash

bam=../data/test.bam
target=chr17
pos=7674420
formats=( fasta fastq )
bindir=../bin
outdir=out/fasta
testdir=ans/fasta

mkdir -p ${outdir}

for format in ${formats[@]}; do
	$bindir/hts-fasta $bam $target $pos A 0 0 $format > ${outdir}/out_${target}-${pos}-A-nodup-nomate.${format}
	$bindir/hts-fasta $bam $target $pos A 1 0 $format > ${outdir}/out_${target}-${pos}-A-dup-nomate.${format}
	$bindir/hts-fasta $bam $target $pos G 0 0 $format > ${outdir}/out_${target}-${pos}-G-nodup-nomate.${format}
	$bindir/hts-fasta $bam $target $pos G 1 0 $format > ${outdir}/out_${target}-${pos}-G-dup-nomate.${format}

	$bindir/hts-fasta $bam $target $pos A 0 1 $format > ${outdir}/out_${target}-${pos}-A-nodup-mate.${format}
	$bindir/hts-fasta $bam $target $pos A 1 1 $format > ${outdir}/out_${target}-${pos}-A-dup-mate.${format}
	$bindir/hts-fasta $bam $target $pos G 0 1 $format > ${outdir}/out_${target}-${pos}-G-nodup-mate.${format}
	$bindir/hts-fasta $bam $target $pos G 1 1 $format > ${outdir}/out_${target}-${pos}-G-dup-mate.${format}
done


echo "Testing diff ...";

dups=( dup nodup )
mates=( mate nomate )

for format in ${formats[@]}; do
	for dup in ${dups[@]}; do
		for mate in ${mates[@]}; do
			diff ${outdir}/out_${target}-${pos}-A-${dup}-${mate}.${format} \
				${testdir}/test_${target}-${pos}-A-${dup}-${mate}.${format}
			echo "${outdir}/out_${target}-${pos}-A-${dup}-${mate}.${format} ... $?"
		done
	done
done
