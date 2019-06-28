#!/bin/bash

bindir=../bin

mkdir -p out/mlat-read-stats

$bindir/hts-mlat-read-stats ans/mlat-read-stats/snvs.tsv ../data/test.bam ../data/tp53_hg38_rc.2bit \
	out/mlat-read-stats/mlat.tsv 7668402 1 0

echo "Testing diff ...";

diff out/mlat-read-stats/mlat.tsv ans/mlat-read-stats/mlat.tsv
echo "out/mlat-read-stats/mlat.tsv ... $?"
