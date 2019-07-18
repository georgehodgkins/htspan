#!/bin/bash

# Merge fastq files (really, any files, other than the line number multiplier) together in integer percentages

set -euo pipefail

if [[ $# -lt 3 ]]; then
	echo "Usage: $0 <dest.fq> <src.fq> <pct_to_merge> [out.fq]"
	exit 1
fi

dest="$1"
src="$2"
pct=$3

if [[ $# -gt 3 ]]; then
	out="$4"
else
	out="$src"
fi

if [[ $pct -gt 100 ]] || [[ $pct -lt 0 ]]; then
	echo "pct_to_merge must be between 0 and 100."
	exit 1
fi

# each fastq record consists of 4 lines
src_lines=$(wc -l < $src)
if [[ $((src_lines % 4)) -ne 0 ]]; then
	echo "Warning: Number of lines in $src is not divisible by 4; is it a FASTQ file?"
fi
src_recs=$((src_lines / 4))

#calculate number of lines to append
recs_app=$(echo "x = $src_recs * ($pct / 100); scale = 0; x / 1" | bc -l)
lines_app=$((recs_app * 4))
echo "Merging $recs_app records ($lines_app lines) from $src and $dest into $out."

head -n -$lines_app $dest > $out
tail -n $lines_app $src >> $out