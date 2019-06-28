ref = ../tp53_hg38.fasta
seed = 3


all: sim.orig.calls.snv sim.ffpe.calls.snv sim.oxog.calls.snv
	

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


sim.orig.calls.snv: sim.orig.bam
	./call.sh $< $(ref)

sim.ffpe.calls.snv: sim.ffpe.bam
	./call.sh $< $(ref)

sim.oxog.calls.snv: sim.oxog.bam
	./call.sh $< $(ref)

diff:
	diff sim.orig.snv sim.orig.calls.snv || true
	diff sim.orig.snv sim.ffpe.calls.snv || true
	diff sim.orig.snv sim.oxog.calls.snv || true

clean:
	rm -rf *.fq *.bam *.bai *.snv

