"""
Quick utility to count how many samples have >= Nx coverage per-variant
"""
import sys
import bisect
import logging
import argparse

import pysam
import truvari

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="dpcnt", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, default="/dev/stdin",
                        help="VCF to annotate (stdin)")
    parser.add_argument("-b", "--bins", type=str, default="0,5,10,15",
                        help="Coverage bins to bisect left the counts (%(default)s)")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output filename (stdout)")
    truvari.setup_logging()
    return parser.parse_args(args)


def edit_header(my_vcf, bins):
    """
    Add INFO for new field to vcf
    """
    header = my_vcf.header.copy()
    desc = 'Description=" Count of samples with >= ' + ",".join([f"{_}x" for _ in bins]) + '">'
    header.add_line('##INFO=<ID=DPCNT,Number=.,Type=Integer,' + desc)
    return header

def add_dpcnt(vcf, n_header=None, bins=None):
    """
    Adds DPCNT to each entry in VCF and yields them
    """
    if bins is None:
        bins = [0, 5, 10, 15]

    if n_header is None:
        n_header = edit_header(vcf, bins)

    for entry in vcf:
        dat = [0] * (len(bins) - 1)
        for sample in entry.samples.values():
            if "DP" in sample and not sample["DP"] is None:
                dp = sample["DP"]
                p = bisect.bisect(bins, dp)
                p = max(0, p - 1)
                dat[p] += 1
        entry.translate(n_header)
        entry.info["DPCNT"] = dat
        yield entry

def dpcnt_main(cmdargs):
    """
    Main method
    """
    args = parse_args(cmdargs)
    try:
        bins = [int(_) for _ in args.bins.split(',')]
    except ValueError:
        logging.error("Cannot parse bins %s", args.bins)

    vcf = pysam.VariantFile(args.input)
    bins.append(sys.maxsize)
    n_header = edit_header(vcf, bins)
    out = pysam.VariantFile(args.output, 'w', header=n_header)
    for entry in add_dpcnt(vcf, n_header, bins):
        out.write(entry)
    out.close()
    logging.info("Finished dpcnt")
