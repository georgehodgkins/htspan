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

if [ $# -ge 4 ]; then
	PASS_ARGS=$4
else
	PASS_ARGS=""
fi

${BIN}/hts-orient-bias identify --stdout -t ${DTYPE} -b ${SIM}/sim.${SRC}.bam -V ${SIM}/sim.${SRC}.calls.snv -o ${SIM}/out/${SRC}.${DTYPE}.${MODEL}.snv --log-file ${SIM}/out/${SRC}.${DTYPE}.${MODEL}.out -M ${MODEL} ${PASS_ARGS}

if [ -e ${SIM}/out/${SRC}.${DTYPE}.${MODEL}.snv ]; then
	echo -e -n "\nFalse positives: "
	echo $(diff ${SIM}/sim.orig.snv ${SIM}/out/${SRC}.${DTYPE}.${MODEL}.snv | grep "^>" | wc -l)

	echo -n "False negatives: "
	echo $(diff ${SIM}/sim.orig.snv ${SIM}/out/${SRC}.${DTYPE}.${MODEL}.snv | grep "^<" | wc -l)
	exit 0
else
	exit 1
fi