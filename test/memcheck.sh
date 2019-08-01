#!/bin/bash

# This file checks for memory leaks on the four separate components of the hts-orient-bias program
# Bayesian quant check may take hours, frequentist quant may take a couple minutes. Ident checks shoud take <30 secs apc.

# check frequentist ident
valgrind --leak-check=summary ../bin/hts-orient-bias identify -t ffpe -M freq -b ../data/sim/sim.ffpe.bam -V ../data/sim/sim.ffpe.calls.vcf -o /dev/null -O vcf --plain --no-stdout
echo

# check Bayesian ident
valgrind --leak-check=summary ../bin/hts-orient-bias identify -t ffpe -M bayes -b ../data/sim/sim.ffpe.bam -V ../data/sim/sim.ffpe.calls.vcf -o /dev/null -O vcf --plain --no-stdout
echo

# check frequentist quant
valgrind --leak-check=summary ../bin/hts-orient-bias quantify -t ffpe -M freq -f ../data/tp53_hg38.fasta -b ../data/sim/sim.ffpe.bam --plain --no-stdout
echo "NOTE: Bayesian quant check takes a /very/ long time. Hope you brought a book or three."

# check Bayesian quant (NB: takes a very long time)
valgrind --leak-check=summary ../bin/hts-orient-bias quantify -t ffpe -M bayes -f ../data/tp53_hg38.fasta -b ../data/sim/sim.ffpe.bam --plain --no-stdout
echo