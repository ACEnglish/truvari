"""
Annotate low complexity region entropy score for variants
Credit: https://jszym.com/blog/dna_protein_complexity/
"""
import math
import logging
import argparse

import pysam

import truvari

def sequence_to_repvec(sequence, N):
    """
    Computes the repetition vector (as seen in Wooton, 1993) from a
    given sequence of a biopolymer with `N` possible residues.

    :param sequence: the nucleotide or protein sequence to generate a repetition vector for.
    :param N: the total number of possible residues in the biopolymer `sequence` belongs to.
    """
    encountered_residues = set()
    repvec = []

    for residue in sequence:
        if residue not in encountered_residues:
            residue_count = sequence.count(residue)

            repvec.append(residue_count)

            encountered_residues.add(residue)

        if len(encountered_residues) == N:
            break

    while len(repvec) < N:
        repvec.append(0)

    return sorted(repvec, reverse=True)

def sequence_entropy(sequence, N=4):
    """
    Computes the Shannon Entropy of a given sequence of a
    biopolymer with `N` possible residues. See (Wooton, 1993)
    for more.

    :param sequence: the nucleotide or protein sequence whose Shannon Entropy is to calculated.
    :param N: the total number of possible residues in the biopolymer `sequence` belongs to.
    """
    repvec = sequence_to_repvec(sequence, N)

    L = len(sequence)

    entropy = sum((-1*(n/L)*math.log((n/L), N) for n in repvec if n != 0))

    return entropy

def edit_header(my_vcf):
    """
    Add INFO for new field to vcf
    """
    header = my_vcf.header.copy()
    header.add_line(('##INFO=<ID=LCR,Number=1,Type=Float,'
                     'Description="Low complexity region entropy score">'))
    return header

def add_lcr(vcf, n_header):
    """
    Adds LCR to each entry in VCF and yields them
    """
    if not n_header:
        n_header = edit_header(vcf)
    for entry in vcf:
        seq = entry.alts[0] if len(entry.alts[0]) > len(entry.ref) else entry.ref
        score = sequence_entropy(seq)
        entry.translate(n_header)
        entry.info["LCR"] = score
        yield entry

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="lcr", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, default="/dev/stdin",
                        help="VCF to annotate (stdin)")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output filename (stdout)")
    truvari.setup_logging()
    return parser.parse_args(args)

def lcr_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)
    vcf = pysam.VariantFile(args.input)
    n_header = edit_header(vcf)
    out = pysam.VariantFile(args.output, 'w', header=n_header)
    for entry in add_lcr(vcf, n_header):
        out.write(entry)
    out.close()
    logging.info("Finished lcr")
