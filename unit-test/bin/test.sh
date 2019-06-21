#!/bin/bash

opts=""

while [ "$1" != "" ]; do
	opts="$opts $1"
	shift
done

./test-bayes-orient-bias-filter $opts
./test-freq-orient-bias-filter $opts
./test-orient-bias-quant $opts
./test-io $opts