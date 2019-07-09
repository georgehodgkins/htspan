#!/bin/bash

DTYPES=("oxog" "ffpe")
TRUE_PHI=".2"
REF_FILE="../data/tp53_hg38.fasta"

declare -a ESTIMATES

echo -e "\tLEGEND:
\t(DTYPE)-c: Control test, testing the undamaged data for DTYPE damage
\tfreq-u: Frequentist test with unknown phi
\tfreq-g: Frequentist test with phi fixed at the true value
\tfreq-e: Frequentist test with phi estimated and fixed at that value\n"

#Header line
echo -e "Dtype\tModel\tFP\tFN"

for DTYPE in "${DTYPES[@]}"; do 
	# controls (test original for DTYPE with the two models)
	RESLINE=$(./ident-test.sh orig ${DTYPE} bayes --no-stdout -v 2)
	echo -e "${DTYPE}-c\tbayes\t${RESLINE% *}\t${RESLINE#* }"

	RESLINE=$(./ident-test.sh orig ${DTYPE} freq --no-stdout -v 2)
	echo -e "${DTYPE}-c\tfreq\t${RESLINE% *}\t${RESLINE#* }"

	# Bayesian
	RESLINE=$(./ident-test.sh ${DTYPE} ${DTYPE} bayes --no-stdout -v 2)
	echo -e "${DTYPE}\tbayes\t${RESLINE% *}\t${RESLINE#* }"

	# Frequentist with unknown phi
	RESLINE=$(./ident-test.sh ${DTYPE} ${DTYPE} freq --no-stdout -v 2)
	echo -e "${DTYPE}\tfreq-u\t${RESLINE% *}\t${RESLINE#* }"

	# Frequentist with phi fixed at ground truth
	RESLINE=$(./ident-test.sh ${DTYPE} ${DTYPE} freq --no-stdout -v 2 -p ${TRUE_PHI} --fixed-phi)
	echo -e "${DTYPE}\tfreq-g\t${RESLINE% *}\t${RESLINE#* }"

	# Frequentist with phi estimated with quantification and fixed at that estimate
	ESTLINE=$(../bin/hts-orient-bias quantify -t ${DTYPE} -f ${REF_FILE} -b ../data/sim/sim.${DTYPE}.bam)
	ESTIMATE=${ESTLINE##* }
	ESTIMATES=("${ESTIMATES[@]}" "${ESTIMATE}")

	RESLINE=$(./ident-test.sh ${DTYPE} ${DTYPE} freq --no-stdout -v 2 -p ${ESTIMATE} --fixed-phi)
	echo -e "${DTYPE}\tfreq-e\t${RESLINE% *}\t${RESLINE#* }"
done

echo -e "\nPhi estimates from quantification:\noxoG: ${ESTIMATES[0]}\nFFPE: ${ESTIMATES[1]}"
exit 0