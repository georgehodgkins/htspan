#!/bin/bash

mkdir -p out/mlat-read-stats

./hts-mlat-read-stats test/mlat-read-stats/snvs.tsv ../data/test.bam ../data/TP53_hg38_rc.2bit \
	out/mlat-read-stats/mlat.tsv 7668402 1 0

echo "Testing diff ...";

diff out/mlat-read-stats/mlat.tsv test/mlat-read-stats/mlat.tsv
echo "out/mlat-read-stats/mlat.tsv ... $?"
