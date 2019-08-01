#!/bin/bash

obias=../../bin/hts-orient-bias
ref=../tp53_hg38.fasta

tags=( 01_01_01 01_01_02 01_01_03 01_01_04 01_01_05 )
#tags=( 01_01_02 01_01_03 01_01_04 01_01_05 )

for tag in "${tags[@]}"; do

	echo $tag

	indir=grid/${tag}
	#tag=${indir##*/}
	outdir=grid/obias/${tag}

	mkdir -p $outdir

	$obias quantify -M freq -t ffpe -f $ref -b $indir/ffpe.bam -J $outdir/ffpe.obquant.json
	phi=$(grep phi $outdir/ffpe.obquant.json | sed -E 's/.*"phi":.([0-9.e+-]+),?/\1/')

	$obias quantify -M bayes -t ffpe -f $ref -b $indir/ffpe.bam --max-reads 20000 -J $outdir/ffpe.bobquant.json
	alpha_phi=$(grep alpha_phi $outdir/ffpe.bobquant.json | sed -E 's/.*"alpha_phi":.([0-9.e+-]+),?/\1/')
	beta_phi=$(grep beta_phi $outdir/ffpe.bobquant.json | sed -E 's/.*"beta_phi":.([0-9.e+-]+),?/\1/')

	# fixed phi
	$obias identify -M freq -t ffpe -b $indir/ffpe.bam -V $indir/ffpe.calls.snv -g 0 -v 0 \
		--phi $phi --fixed-phi -o $outdir/ffpe.fixed.snv

	# unknown phi
	$obias identify -M freq -t ffpe -b $indir/ffpe.bam -V $indir/ffpe.calls.snv -g 0 -v 0 \
		--phi $phi -o $outdir/ffpe.unknown.snv

	# variable phi
	$obias identify -M bayes -t ffpe -b $indir/ffpe.bam -V $indir/ffpe.calls.snv -g 0 -v 0 \
		--alpha $alpha_phi --beta $beta_phi -o $outdir/ffpe.variable.snv

done
