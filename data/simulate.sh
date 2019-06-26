
ref=tp53_hg38.fasta
out=sim

r1=$out.r1.fq
r2=$out.r2.fq

# 2*100 bp reads * 10000 reads / 19149 bp seq = 104x coverage
wgsim -h -1 100 -2 100 -S 0 -e 0.05 -N 10000 -R 0 \
	$ref $r1 $r2 | 
	cut -f 1-4 > $out.snv

# Align reads

if [[ ! -f ${ref}.bwt ]]; then
	bwa index ${ref}
fi

bwa mem $ref $r1 $r2 | 
	samtools view -b > $out.raw.bam
samtools sort $out.raw.bam -o $out.bam
rm $out.raw.bam

