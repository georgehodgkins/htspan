#!/bin/bash

mkdir -p out/mlat-filter

./hts-mlat-filter test/mlat-read-stats/snvs.tsv ../data/test.bam ../data/TP53_hg38_rc.2bit \
	out/mlat-filter/pass.vtr 0.5 7668402

echo "Testing diff ...";

diff out/mlat-filter/pass.vtr test/mlat-filter/pass.vtr
echo "out/mlat-filter/pass.vtr ... $?"
