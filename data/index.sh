#!/bin/bash

ref=tp53_hg38.fasta

bwa index $ref
samtools faidx $ref
gatk CreateSequenceDictionary -R $ref -O ${ref%.*}
