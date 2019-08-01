#!/bin/bash

# replace e with 10^ so bc can interpret scientific notation
sci_to_dec () {
	echo $(sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<< "$1")
}

set -euxo pipefail

SIM="../data/sim"
BIN="../bin"
REF="../data/tp53_hg38.fasta"

# maximum estimate for a no damage case (i.e. checking oxoG damaged or original for FFPE) to pass
ND_MAX=".01"
# minimum and maximum estimates for a damage case (i.e. checking oxoG damaged for oxoG)
DM_MIN=".18"
DM_MAX=".22"

# turn off command reporting after setting initial parameters
set +x

# build program
make -C .. bin/hts-orient-bias

# build sim data
make -C $SIM

# obtain estimates
ORIG_OXOG="$(sci_to_dec $($BIN/hts-orient-bias quantify -M freq -t oxog -f $REF -b $SIM/sim.orig.bam --plain))"
ORIG_FFPE="$(sci_to_dec $($BIN/hts-orient-bias quantify -M freq -t ffpe -f $REF -b $SIM/sim.orig.bam --plain))"

FFPE_OXOG="$(sci_to_dec $( $BIN/hts-orient-bias quantify -M freq -t oxog -f $REF -b $SIM/sim.ffpe.bam --plain))"
FFPE_FFPE="$(sci_to_dec $( $BIN/hts-orient-bias quantify -M freq -t ffpe -f $REF -b $SIM/sim.ffpe.bam --plain))"

OXOG_OXOG="$(sci_to_dec $( $BIN/hts-orient-bias quantify -M freq -t oxog -f $REF -b $SIM/sim.oxog.bam --plain))"
OXOG_FFPE="$(sci_to_dec $( $BIN/hts-orient-bias quantify -M freq -t ffpe -f $REF -b $SIM/sim.oxog.bam --plain))"

FAIL=0

# report on estimates
if [ $(echo "$ORIG_OXOG < $ND_MAX" | bc) -eq 0 ]; then
	echo "Estimate of oxoG damage on undamaged BAM was too high. Got: $ORIG_OXOG"
	FAIL=1
fi

if [ $(echo "$ORIG_FFPE< $ND_MAX" | bc) -eq 0 ]; then
	echo "Estimate of FFPE damage on undamaged BAM was too high. Got: $ORIG_FFPE"
	FAIL=1
fi

if [ $(echo "$FFPE_OXOG < $ND_MAX" | bc) -eq 0 ]; then
	echo "Estimate of oxoG damage on FFPE-damaged BAM was too high. Got: $FFPE_OXOG"
	FAIL=1
fi

if [ $(echo "$FFPE_FFPE < $DM_MAX" | bc) -eq 0 ]; then
	echo "Estimate of FFPE damage on FFPE-damaged BAM was too high. Got: $FFPE_FFPE"
	FAIL=1
elif [ $(echo "$FFPE_FFPE > $DM_MIN" | bc) -eq 0 ]; then
	echo "Estimate of FFPE damage on FFPE-damaged BAM was too low. Got: $FFPE_FFPE"
	FAIL=1
fi

if [ $(echo "$OXOG_OXOG < $DM_MAX" | bc) -eq 0 ]; then
	echo "Estimate of oxoG damage on oxoG-damaged BAM was too high. Got: $OXOG_OXOG"
	FAIL=1
elif [ $(echo "$OXOG_OXOG > $DM_MIN" | bc) -eq 0 ]; then
	echo "Estimate of oxoG damage on oxoG-damaged BAM was too low. Got: $OXOG_OXOG"
	FAIL=1
fi

if [ $(echo "$OXOG_FFPE < $ND_MAX" | bc) -eq 0 ]; then
	echo "Estimate of FFPE damage on oxoG-damaged BAM was too high. Got: $OXOG_FFPE"
	FAIL=1
fi

if [ $FAIL -eq 0 ]; then
	echo "All tests passed."
fi