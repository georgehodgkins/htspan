#!/bin/bash

BIN="../bin"
SIM="../data/sim"

if [ $# -ge 3 ]; then
	SRC=$1
	DTYPE=$2
	MODEL=$3
else
	echo "Arguments: <src=orig|ffpe|oxog> <dtype=ffpe|oxog> <model=bayes|freq> [<args_to_pass>]"
	exit 1
fi

PASS_ARGS=""


if [ $# -ge 4 ]; then
	while [ "$4" != "" ]; do
		PASS_ARGS="$PASS_ARGS $4"
		shift
	done
fi

rm -f ${SIM}/out/${SRC}.${DTYPE}.${MODEL}.*

${BIN}/hts-orient-bias identify --stdout -t ${DTYPE} -b ${SIM}/sim.${SRC}.bam -V ${SIM}/sim.${SRC}.calls.snv -o out/${SRC}.${DTYPE}.${MODEL}.snv --log-file out/${SRC}.${DTYPE}.${MODEL}.out -M ${MODEL} ${PASS_ARGS}

if [ -e out/${SRC}.${DTYPE}.${MODEL}.snv ]; then
	echo -n $(diff ${SIM}/sim.orig.snv <(cut -f -4 out/${SRC}.${DTYPE}.${MODEL}.snv) | grep "^>" | wc -l) #false positives
	echo -n " "
	echo $(diff ${SIM}/sim.orig.snv <(cut -f -4 out/${SRC}.${DTYPE}.${MODEL}.snv) | grep "^<" | wc -l) #false negatives
	exit 0
else
	exit 1
fi