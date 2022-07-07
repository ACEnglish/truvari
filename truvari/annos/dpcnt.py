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
    parser.add_argument("--no-ad", action="store_false",
                        help="Skip adding ADCNT bins")
    parser.add_argument("-p", "--present", action="store_true", default=False,
                        help="Only count sites with present (non ./.) genotypes")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output filename (stdout)")
    truvari.setup_logging()
    return parser.parse_args(args)


def edit_header(my_vcf, bins, add_ad=False):
    """
    Add INFO for new field to vcf
    """
    header = my_vcf.header.copy()
    desc = 'Description="Count of samples with >= ' + ",".join([f"{_}x" for _ in bins][:-1]) + '">'
    header.add_line('##INFO=<ID=DPCNT,Number=.,Type=Integer,' + desc)
    if add_ad:
        header.add_line('##INFO=<ID=ADCNT,Number=.,Type=Integer,' + desc)
    return header

def add_dpcnt(vcf, n_header=None, bins=None, add_ad=False, present=False):
    """
    Adds DPCNT to each entry in VCF and yields them
    """
    if bins is None:
        bins = [0, 5, 10, 15, sys.maxsize]

    if n_header is None:
        n_header = edit_header(vcf, bins, add_ad)

    for entry in vcf:
        dat = [0] * (len(bins) - 1)
        dat_ad = [0] * (len(bins) - 1)
        for sample in entry.samples.values():
            if present and sample["GT"] == (None, None):
                continue
            if "DP" in sample and not sample["DP"] is None:
                dp = sample["DP"]
                p = bisect.bisect(bins, dp)
                p = max(0, p - 1)
                dat[p] += 1
            if add_ad and "AD" in sample and len(sample["AD"]) == 2:
                ad = sample["AD"][1]
                if isinstance(ad, int):
                    p = bisect.bisect(bins, ad)
                    p = max(0, p - 1)
                    dat_ad[p] += 1

        entry.translate(n_header)
        entry.info["DPCNT"] = dat
        if add_ad:
            entry.info["ADCNT"] = dat_ad
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
    n_header = edit_header(vcf, bins, args.no_ad)
    out = pysam.VariantFile(args.output, 'w', header=n_header)
    for entry in add_dpcnt(vcf, n_header, bins, args.no_ad, args.present):
        out.write(entry)
    out.close()
    logging.info("Finished dpcnt")
