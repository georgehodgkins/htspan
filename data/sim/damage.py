#!/usr/bin/env python3

import sys
import argparse
import random

class fastq_record:
    """A FASTQ record."""
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual

def read_fastq_record(f):
    """Read a FASTQ record from file f or return None on failure."""
    name = f.readline().rstrip()
    if not name: return None

    seq = f.readline().rstrip()
    if not seq: return None

    marker = f.readline().rstrip()
    if marker != "+":
        print("ERROR: FASTQ file is malformed", stderr)
        return None

    qual = f.readline().rstrip()
    if not qual: return None

    return fastq_record(name, seq, qual)


def write_fastq_record(r, f):
    """Write a FASTQ record r to file f."""
    f.write('\n'.join([r.name, r.seq, "+", r.qual, ""]))


def damage_seq(xs, ref, alt, phi):
    """
    Randomly mutate ref nucleotides in sequence xs to alt nucleotides
    with probability phi.
    """
    ys = ""
    for n in xs:
        if n.upper() == ref:
            if random.random() < phi:
                ys += alt
                continue
        ys += n
    return ys

def prepend_ext(x, tag):
    """Insert file extension tag after first period."""
    i = x.find('.')
    return x[:i] + "." + tag + x[i:]


if __name__ == "__main__":

    pr = argparse.ArgumentParser(description = "Introduce damage into Illumina sequencing data")
    pr.add_argument("r1", help="fastq file containing R1 reads")
    pr.add_argument("r2", help="fastq file containing R2 reads")
    pr.add_argument("--out1", help="output file name for R1 reads")
    pr.add_argument("--out2", help="output file name for R2 reads")
    pr.add_argument("--phi", type=float, default=0.1, help="damage rate")
    pr.add_argument("--seed", type=int, default=0, help="random seed")
    pr.add_argument("--type", default="FFPE", help="damage type [FFPE or oxoG]")

    args = pr.parse_args()

    phi = args.phi

    # set random seed to fixed or arbitrary value
    if args.seed != 0:
        random.seed(args.seed)
    else:
        random.seed()

    # configure damage type
    dtype = args.type.lower()
    if dtype == "oxog":
        nt_ref_r1 = 'G'; nt_alt_r1 = 'T'
        nt_ref_r2 = 'C'; nt_alt_r2 = 'A'
    elif dtype == "ffpe":
        nt_ref_r1 = 'C'; nt_alt_r1 = 'T'
        nt_ref_r2 = 'G'; nt_alt_r2 = 'A'
    else:
        print("ERROR: Damage type is not supported", stderr)

    # set output file names
    if not args.out1:
        out1 = prepend_ext(args.r1, dtype)
    if not args.out2:
        out2 = prepend_ext(args.r2, dtype)

    # read input sequences and write damaged sequences
    with open(args.r1) as r1f, open(args.r2) as r2f, \
         open(out1, 'w') as out1f, open(out2, 'w') as out2f:
        
        while True:
            r1 = read_fastq_record(r1f)
            r2 = read_fastq_record(r2f)
            if not r1 or not r2: break

            r1_new_seq = damage_seq(r1.seq, nt_ref_r1, nt_alt_r1, phi)
            r2_new_seq = damage_seq(r2.seq, nt_ref_r2, nt_alt_r2, phi)

            s1 = fastq_record(r1.name, r1_new_seq, r1.qual)
            s2 = fastq_record(r2.name, r2_new_seq, r2.qual)
            
            write_fastq_record(s1, out1f)
            write_fastq_record(s2, out2f)

