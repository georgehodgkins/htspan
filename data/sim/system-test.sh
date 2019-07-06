#!/bin/bash

BIN="../../bin"

#cd ../.. && make clean && make && cd data/sim

# Control (testing undamaged BAM for damage)
${BIN}/hts-orient-bias identify -t oxog -M freq -b sim.orig.bam -V sim.oxog.calls.snv -o out/orig.oxog.freq.snv --result-file out/orig.oxog.freq.vals
${BIN}/hts-orient-bias identify -t oxog -M bayes -b sim.orig.bam -V sim.oxog.calls.snv -o out/orig.oxog.bayes.snv --result-file out/orig.oxog.bayes.vals
${BIN}/hts-orient-bias identify -t ffpe -M freq -b sim.orig.bam -V sim.ffpe.calls.snv -o out/orig.ffpe.freq.snv --result-file out/orig.ffpe.freq.vals
${BIN}/hts-orient-bias identify -t ffpe -M bayes -b sim.orig.bam -V sim.ffpe.calls.snv -o out/orig.ffpe.bayes.snv --result-file out/orig.ffpe.bayes.vals

# oxoG damaged reads
${BIN}/hts-orient-bias identify -t oxog -M freq -b sim.oxog.bam -V sim.oxog.calls.snv -o out/oxog.oxog.freq.snv --result-file out/oxog.oxog.freq.vals
${BIN}/hts-orient-bias identify -t oxog -M bayes -b sim.oxog.bam -V sim.oxog.calls.snv -o out/oxog.oxog.bayes.snv --result-file out/oxog.oxog.bayes.vals
${BIN}/hts-orient-bias identify -t ffpe -M freq -b sim.oxog.bam -V sim.ffpe.calls.snv -o out/oxog.ffpe.freq.snv --result-file out/oxog.ffpe.freq.vals
${BIN}/hts-orient-bias identify -t ffpe -M bayes -b sim.oxog.bam -V sim.ffpe.calls.snv -o out/oxog.ffpe.bayes.snv --result-file out/oxog.ffpe.bayes.vals

# FFPE damaged reads
${BIN}/hts-orient-bias identify -t oxog -M freq -b sim.ffpe.bam -V sim.oxog.calls.snv -o out/ffpe.oxog.freq.snv --result-file out/ffpe.oxog.freq.vals
${BIN}/hts-orient-bias identify -t oxog -M bayes -b sim.ffpe.bam -V sim.oxog.calls.snv -o out/ffpe.oxog.bayes.snv --result-file out/ffpe.oxog.bayes.vals
${BIN}/hts-orient-bias identify -t ffpe -M freq -b sim.ffpe.bam -V sim.ffpe.calls.snv -o out/ffpe.ffpe.freq.snv --result-file out/ffpe.ffpe.freq.vals
${BIN}/hts-orient-bias identify -t ffpe -M bayes -b sim.ffpe.bam -V sim.ffpe.calls.snv -o out/ffpe.ffpe.bayes.snv --result-file out/ffpe.ffpe.bayes.vals