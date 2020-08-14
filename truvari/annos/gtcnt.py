"""
Annotates GTCounts of alleles
"""
import argparse

import pysam
import truvari

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="gtcnt", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, default="/dev/stdin",
                        help="VCF to annotate (stdin)")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output filename (stdout)")
    return parser.parse_args(args)


def edit_header(my_vcf):
    """
    Add INFO for new field to vcf
    """
    header = my_vcf.header.copy()
    header.add_line(('##INFO=<ID=GTCNT,Number=.,Type=Integer,'
                     'Description="Counts of genotypes for the allele (UNK, REF, HET, HOM)">'))
    return header

def add_gtcnt(vcf, out, n_header=None):
    """
    Adds GTCNT to each entry in VCF and yields them
    """
    if n_header is None:
        n_header = edit_header(vcf)
    for entry in vcf:
        cnt = [0, 0, 0, 0]
        #cnt = {"UNK": 0, "REF": 0, "HET": 0, "HOM": 0}
        for sample in entry.samples:
            gt = entry.samples[sample]["GT"]
            if None in gt or len(gt) != 2:
                cnt[0] += 1
            elif gt[0] == gt[1] and gt[0] == 0:
                cnt[1] += 1
            elif gt[0] == gt[1]:
                cnt[3] += 1
            elif gt[0] != gt[1]:
                cnt[2] += 1
            else:
                cnt[0] += 1
        try:
            nentry = truvari.copy_entry(entry, n_header)
        except TypeError:
            yield entry
            continue
        nentry.info["GTCNT"] = cnt 
        yield nentry

def gtcnt_main(cmdargs):
    """
    Main method
    """
    args = parse_args(cmdargs)
    vcf = pysam.VariantFile(args.input)
    n_header = edit_header(vcf)
    out = pysam.VariantFile(args.output, 'w', header=n_header)
    for entry in add_gtcnt(vcf, out, n_header):
        out.write(entry)
    out.close()
