bwa index TP53_hg38.fasta
bwa mem TP53_hg38.fasta reads_tp53.fastq | samtools view -b > reads_tp53.raw.bam
samtools sort reads_tp53.raw.bam -o reads_tp53.bam
rm reads_tp53.raw.bam

