"""
Annotates GC content of SVs
"""
import re
import sys
import argparse

import pysam
import pyfaidx
import truvari

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="gcpct", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, default="/dev/stdin",
                        help="VCF to annotate (stdin)")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output filename (stdout)")
    parser.add_argument("-r", "--reference", type=str, required=True,
                        help="Reference fasta")
    return parser.parse_args(args)


def edit_header(my_vcf):
    """
    Add INFO for new fields to vcf
    #Probably want to put in the PG whatever, too
    """
    # Update header
    # Edit Header
    header = my_vcf.header.copy()
    header.add_line(('##INFO=<ID=GCPCT,Number=1,Type=Integer,'
                     'Description="GC Percent of the reference call range or alt sequence (whichever is longer)">'))
    return header

def add_gcpct(vcf, ref, out, n_header=None):
    """
    Adds GCPCT to each entry in VCF and yields them
    """
    if not n_header:
        n_header = edit_header(vcf)

    for entry in vcf:
        start, end = truvari.entry_boundaries(entry)
        span = abs(end - start)
        try:
            seq = ref.get_seq(entry.chrom, start, end).seq if span >= len(entry.alts[0]) else str(entry.alts[0])
            gcpct = int((sum(1 for m in re.finditer("[GC]", seq)) / len(seq)) * 100)
            entry = truvari.copy_entry(entry, n_header)
            entry.info["GCPCT"] = gcpct
        except Exception:
            # Silent failures aren't the best
            pass
        yield entry

def gcpct_main(cmdargs):
    """
    Main method
    """
    args = parse_args(cmdargs)
    ref = pyfaidx.Fasta(args.reference)
    vcf = pysam.VariantFile(args.input)
    n_header = edit_header(vcf)
    out = pysam.VariantFile(args.output, 'w', header=n_header)
    for entry in add_gcpct(vcf, ref, out, n_header):
        out.write(entry)
    out.close()
