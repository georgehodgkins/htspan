#!/bin/bash

vaf=../../bin/hts-vaf
ref=../tp53_hg38.fasta


#dset=grid1
#groups=( 01_01_01 01_01_02 01_01_03 01_01_04 01_01_05 05_01_01 05_01_02 05_01_03 05_01_04 05_01_05 )

#dset=grid2
#groups=( 01_01_01 01_01_02 01_01_03 01_01_04 05_01_01 05_01_02 05_01_03 05_01_04 )

dset=grid3
groups=( 01_01_01 01_01_02 01_01_03 01_01_04 05_01_01 05_01_02 05_01_03 05_01_04 )

#dset=grid4
#groups=( 01_01_01 01_01_02 01_01_03 01_01_04 05_01_01 05_01_02 05_01_03 05_01_04 )


for group in "${groups[@]}"; do

	echo $group

	indir="${dset}/input/${group}"
	outdir="${dset}/vaf/${group}"

	mkdir -p $outdir

	$vaf C T $indir/ffpe.bam $indir/ffpe.calls.snv $outdir/ffpe.vaf.snv

done

