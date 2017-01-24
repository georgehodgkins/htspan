#!/bin/bash

bam=../data/test.bam
target=chr17
pos=7674420
formats=( fasta fastq )

mkdir -p out

for format in ${formats[@]}; do
	./hts-fasta $bam $target $pos A 0 0 $format > out/out_${target}-${pos}-A-nodup-nomate.${format}
	./hts-fasta $bam $target $pos A 1 0 $format > out/out_${target}-${pos}-A-dup-nomate.${format}
	./hts-fasta $bam $target $pos G 0 0 $format > out/out_${target}-${pos}-G-nodup-nomate.${format}
	./hts-fasta $bam $target $pos G 1 0 $format > out/out_${target}-${pos}-G-dup-nomate.${format}

	./hts-fasta $bam $target $pos A 0 1 $format > out/out_${target}-${pos}-A-nodup-mate.${format}
	./hts-fasta $bam $target $pos A 1 1 $format > out/out_${target}-${pos}-A-dup-mate.${format}
	#./hts-fasta $bam $target $pos G 0 1 $format > out/out_${target}-${pos}-G-nodup-mate.${format}
	#./hts-fasta $bam $target $pos G 1 1 $format > out/out_${target}-${pos}-G-dup-mate.${format}
done


echo "Testing diff ...";

dups=( dup nodup )
mates=( mate nomate )

for format in ${formats[@]}; do
	for dup in ${dups[@]}; do
		for mate in ${mates[@]}; do
			diff \
				out/out_${target}-${pos}-A-${dup}-${mate}.${format} \
				test/test_${target}-${pos}-A-${dup}-${mate}.${format}
			echo "out/out_${target}-${pos}-A-${dup}-${mate}.${format} ... $?"
		done
	done
done
