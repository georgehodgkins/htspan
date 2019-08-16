#!/bin/bash

ref=../tp53_hg38.fasta

dset=grid1

#groups=( 01_01_01 01_01_02 01_01_03 01_01_04 01_01_05 )
groups=( 05_01_01 05_01_02 05_01_03 05_01_04 05_01_05 )

for group in "${groups[@]}"; do

	echo $group

	indir="${dset}/input/${group}"
	outdir="${dset}/gatk/${group}"

	mkdir -p $outdir

	gatk CollectSequencingArtifactMetrics \
		-I $indir/ffpe.bam -R $ref -O $outdir/ffpe

done

