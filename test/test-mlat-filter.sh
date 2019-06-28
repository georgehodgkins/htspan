#!/bin/bash

bindir=../bin

mkdir -p out/mlat-filter

$bindir/hts-mlat-filter ans/mlat-read-stats/snvs.tsv ../data/test.bam ../data/tp53_hg38_rc.2bit \
	out/mlat-filter/pass.vtr 0.5 7668402

echo "Testing diff ...";

diff out/mlat-filter/pass.vtr ans/mlat-filter/pass.vtr
echo "out/mlat-filter/pass.vtr ... $?"
