#!/usr/bin/env python3


if __name__ == "__main__":

    import argparse
    import random
    import subprocess

    pr = argparse.ArgumentParser("Cancer short read simulator")
    pr.add_argument("-T", "--theta", type=float, help="rate of genuine mutations", default=0.02)
    pr.add_argument("-e", "--base_error", type=float, help="rate of base error", default=0.005)
    pr.add_argument("-p", "--purity", type=float, help="cancer purity", default=0.8)
    pr.add_argument("-s", "--seed", type=int, help="seed for random number generator", default=0)
    #pr.add_argument("-q", "--homozygous", help="homozygous mutation probability", default=0.1)
    pr.add_argument("-b", "--wgsim", help="path to wgsim", default="wgsim")
    pr.add_argument("-l", "--read_length", type=int, help="read length", default=100)
    pr.add_argument("-i", "--insert_size", type=int, help="insert size", default=500)
    pr.add_argument("-n", "--nreads", type=int, help="number of read pairs", default=10000)
    pr.add_argument("reference", help="reference sequence")
    pr.add_argument("prefix", help="output file name prefix")

    args = pr.parse_args()

    wgsim = args.wgsim
    reference = args.reference
    prefix = args.prefix

    seed = args.seed
    theta = args.theta
    berror = args.base_error
    purity = args.purity
    #hprob = args.homozygous
    rlen = args.read_length
    isize = args.insert_size
    nreads = args.nreads
    prefix = args.prefix

    # set random seed to fixed or arbitrary value
    if seed != 0:
        random.seed(seed)
    else:
        random.seed()


    cmd_base = [
        wgsim,
        "-e", str(berror),
        "-d", str(isize),
        "-s", str(round(isize * 0.1)),
        "-1", str(rlen), "-2", str(rlen),
        "-R", "0",
        "-S", str(seed),
        "-h"
    ]

    normal_h1_r1 = prefix + ".normal.h1.r1.fq"
    normal_h1_r2 = prefix + ".normal.h1.r2.fq"
    normal_h1_nreads = round(nreads * (1 - purity) * 0.5)
    cmd_normal_h1 = cmd_base + ["-r", "0", "-N", str(normal_h1_nreads)] + [reference, normal_h1_r1, normal_h1_r2]

    normal_h2_r1 = prefix + ".normal.h2.r1.fq"
    normal_h2_r2 = prefix + ".normal.h2.r2.fq"
    normal_h2_nreads = round(nreads * (1 - purity) * 0.5)
    cmd_normal_h2 = cmd_base + ["-r", "0", "-N", str(normal_h2_nreads)] + [reference, normal_h2_r1, normal_h2_r2]

    # first cancer haplotype contains no somatic mutations
    cancer_h1_r1 = prefix + ".cancer.h1.r1.fq"
    cancer_h1_r2 = prefix + ".cancer.h1.r2.fq"
    cancer_h1_nreads = round(nreads * purity * 0.5)
    cmd_cancer_h1 = cmd_base + ["-r", "0", "-N", str(cancer_h1_nreads)] + [reference, cancer_h1_r1, cancer_h1_r2]

    # second cancer haplotype contains somatic mutations
    cancer_h2_r1 = prefix + ".cancer.h2.r1.fq"
    cancer_h2_r2 = prefix + ".cancer.h2.r2.fq"
    cancer_h2_nreads = round(nreads * purity * 0.5)
    cmd_cancer_h2 = cmd_base + ["-r", str(theta), "-N", str(cancer_h2_nreads)] + [reference, cancer_h2_r1, cancer_h2_r2]

    subprocess.run(cmd_normal_h1)
    subprocess.run(cmd_normal_h2)
    subprocess.run(cmd_cancer_h1)

    snv = prefix + ".snv"

    ps = subprocess.Popen(cmd_cancer_h2, stdout=subprocess.PIPE)
    with open(snv, "wb") as outf_snv:
        outf_snv.write(b'\t'.join([b"chrom", b"pos", b"ref", b"alt"]) + b'\n')
        for line in ps.stdout:
            # keep only the first four columns
            x = b'\t'.join(line.rstrip().split(b'\t')[:4]) + b'\n'
            outf_snv.write(x)

    r1 = prefix + ".r1.fq"
    r2 = prefix + ".r2.fq"

    with open(r1, "w") as outf_r1:
        subprocess.run(["cat", normal_h1_r1, normal_h2_r1, cancer_h1_r1, cancer_h2_r1], stdout=outf_r1)

    with open(r2, "w") as outf_r2:
        subprocess.run(["cat", normal_h1_r2, normal_h2_r2, cancer_h1_r2, cancer_h2_r2], stdout=outf_r2)

    subprocess.run(["rm", normal_h1_r1, normal_h2_r1, cancer_h1_r1, cancer_h2_r1])
    subprocess.run(["rm", normal_h1_r2, normal_h2_r2, cancer_h1_r2, cancer_h2_r2])
    

