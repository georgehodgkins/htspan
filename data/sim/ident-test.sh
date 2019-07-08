#!/bin/bash

BIN="../../bin"

if [ $# -ge 3 ]; then
	SRC=$1
	DTYPE=$2
	MODEL=$3
else
	echo "Arguments: <src=orig|ffpe|oxog> <dtype=orig|ffpe|oxog> <model=bayes|freq> [<args_to_pass>]"
	exit 1
fi

if [ $# -ge 4 ]; then
	PASS_ARGS=$4
else
	PASS_ARGS=""
fi

../../bin/hts-orient-bias identify --stdout -t ${DTYPE} -b sim.${SRC}.bam -V sim.${DTYPE}.calls.snv -o out/${SRC}.${DTYPE}.${MODEL}.snv --log-file out/${SRC}.${DTYPE}.${MODEL}.out -M ${MODEL} ${PASS_ARGS}