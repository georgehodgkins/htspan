#!/bin/bash

set -euo pipefail

if [[ $# -ge 3 || $# -ge 1 && $1 == "-h" ]]; then
	echo "usage $0 [in.vcf | in.vcf.gz] [out.snv]" >&2
	exit 1
fi

if [[ $# -ge 1 ]]; then
	infile=$1
else
	infile=/dev/stdin
fi

if [[ $# -ge 2 ]]; then
	outfile=$2
else
	outfile=/dev/stdout
fi

fname=${infile##*/}
fext=${fname##*.}
fstem=${fname%%.*}

outfile=${outfile:-${fstem}.snv}

if [[ $fext == "gz" ]]; then
	cat=zcat
else
	cat=cat
fi

printf 'chrom\tpos\tref\talt\n' > $outfile
$cat $infile | grep -v '^#' | cut -f 1,2,4,5 >> $outfile

