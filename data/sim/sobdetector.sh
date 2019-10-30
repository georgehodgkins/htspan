#!/bin/bash

obias=../../bin/hts-orient-bias
ref=../tp53_hg38.fasta
jar=~/java/SOBDetector_v0.1.jar

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
    outdir="${dset}/sobdetector/${group}"

    mkdir -p $outdir

    java -jar $jar --input-type Table --input-variants $indir/ffpe.calls.snv --input-bam $indir/ffpe.bam --output-variants $outdir/ffpe.sobdetector.snv
    # SOBDetector inserts an extra TAB character after first 4 columns
    sed -i 's/\t\t/\t/g' $outdir/ffpe.sobdetector.snv

done

