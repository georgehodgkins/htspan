#!/bin/bash

mkdir -p out/mlat-filter

./hts-mlat-filter test/mlat-filter/snvs.tsv ../data/test.bam ../data/TP53_hg38_rc.2bit \
	out/mlat-filter/mlat.tsv 7668402 1 0

echo "Testing diff ...";

diff out/mlat-filter/mlat.tsv test/mlat-filter/mlat.tsv
echo "out/mlat-filter/mlat.tsv ... $?"
