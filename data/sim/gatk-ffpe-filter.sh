#!/bin/bash

set -euo pipefail

ref=../tp53_hg38.fasta


#dset=grid1
#groups=( 01_01_01 01_01_02 01_01_03 01_01_04 01_01_05 )
#groups=( 05_01_01 05_01_02 05_01_03 05_01_04 05_01_05 )

#dset=grid2
#groups=( 01_01_01 01_01_02 01_01_03 01_01_04 05_01_01 05_01_02 05_01_03 05_01_04 )

dset=grid3
groups=( 01_01_01 01_01_02 01_01_03 01_01_04 05_01_01 05_01_02 05_01_03 05_01_04 )

#dset=grid4
#groups=( 01_01_01 01_01_02 01_01_03 01_01_04 )
#groups=( 05_01_01 05_01_02 05_01_03 05_01_04 )


for group in "${groups[@]}"; do

	echo $group

	indir="${dset}/input/${group}"
	outdir="${dset}/gatk/${group}"

	mkdir -p $outdir

	gatk CollectSequencingArtifactMetrics \
		-R $ref \
		-I $indir/ffpe.bam \
		-O $outdir/ffpe

	gatk FilterByOrientationBias \
		-V $indir/ffpe.calls.vcf --artifact-modes 'C/T' \
		-P $outdir/ffpe.pre_adapter_detail_metrics \
		-O $outdir/ffpe.gatk.vcf

	./gatk-ob-vcf2snv.R $outdir/ffpe.gatk.vcf -o $outdir/ffpe.gatk.snv

done

