#!/bin/bash

set -euo pipefail

ref=../tp53_hg38.fasta

outdir=grid
home="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
bin=$home

mkdir -p $outdir

if [[ ! -f $ref ]]; then
	echo "Error: $ref does not exist." >&2
	exit 1
fi

berror=0.001
nreads=50000
purities=( 0.90 0.75 0.50 0.25 0.10 )
thetas=( 0.05 )
phis=( 0.4 0.2 0.1 0.05 0.025 )
damages=( ffpe oxog )

seed=1

i=1
for purity in "${purities[@]}"; do
	j=1
	for theta in "${thetas[@]}"; do
		k=1
		for phi in "${phis[@]}"; do
			I=$(printf '%02d' $i)
			J=$(printf '%02d' $j)
			K=$(printf '%02d' $k)

			prefix=$outdir/${I}_${J}_${K}
			mkdir -p $prefix

			params="{ \"purity\": $purity, \"theta\": $theta, \"phi\": $phi }"
			echo $params | tee $prefix/params.json >&2

			# simulate reads
			$bin/csrsim.py -s $seed -e $berror -p $purity -T $theta -n $nreads $ref $prefix/orig

			# align reads
			$bin/align.sh $prefix/orig $ref

			# call snvs
			$bin/freebayes-call.sh $prefix/orig.bam $ref

			for damage in "${damages[@]}"; do
				# damage reads
				$bin/damage.py -s $seed -t $damage -P $phi $prefix/orig.r1.fq $prefix/orig.r2.fq \
					-1 $prefix/${damage}.r1.fq -2 $prefix/${damage}.r2.fq

				# align reads
				$bin/align.sh $prefix/$damage $ref

				# call snvs
				$bin/freebayes-call.sh $prefix/$damage.bam $ref
			done
			((k++))
		done
		((j++))
	done
	((i++))
done

