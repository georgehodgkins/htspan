ref = ../tp53_hg38.fasta
seed = 1


all: sim.orig.bam sim.ffpe.bam sim.oxog.bam
	

sim.orig.snv:
	./simulate.sh sim.orig $(ref) $(seed)

sim.orig.r1.fq: sim.orig.snv
	

sim.orig.r2.fq: sim.orig.snv
	

sim.orig.bam: sim.orig.r1.fq sim.orig.r2.fq
	./align.sh sim.orig $(ref)

sim.ffpe.r1.fq: sim.orig.r1.fq sim.orig.r2.fq
	./damage.py --type ffpe sim.orig.r1.fq sim.orig.r2.fq --out1 sim.ffpe.r1.fq --out2 sim.ffpe.r2.fq --seed $(seed)

sim.ffpe.r2.fq: sim.ffpe.r1.fq
	

sim.ffpe.bam: sim.ffpe.r1.fq sim.ffpe.r2.fq
	./align.sh sim.ffpe $(ref)

sim.oxog.r1.fq: sim.orig.r1.fq sim.orig.r2.fq
	./damage.py --type oxog sim.orig.r1.fq sim.orig.r2.fq --out1 sim.oxog.r1.fq --out2 sim.oxog.r2.fq --seed $(seed)

sim.oxog.r2.fq: sim.oxog.r1.fq
	

sim.oxog.bam: sim.oxog.r1.fq sim.oxog.r2.fq
	./align.sh sim.oxog $(ref)

clean:
	rm -rf *.bam *.fq *.snv

