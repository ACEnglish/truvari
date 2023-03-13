"""
Count variants per-region in vcf
"""
import os
import argparse

import pysam
import pandas as pd

import truvari

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="stratify", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("regions", metavar="BED",
                        help="Regions to process")
    parser.add_argument("in_vcf", metavar="VCF",
                        help="Truvari bench result directory or a single VCF")
    parser.add_argument("-o", "--output", metavar="OUT", default="/dev/stdout",
                        help="Output bed-like file")
    parser.add_argument("--header", action="store_true",
                        help="Input regions have header to preserve in output")
    parser.add_argument("-w", "--within", action="store_true",
                        help="Only count variants contained completely within region boundaries")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug, show_version=True)
    return args

def count_entries(vcf, regions, within):
    """
    Count the number of variants per bed region a VCF
    regions is a list of lists with sub-lists having 0:3 of chrom, start, end
    Returns a list of the counts in the same order as the regions
    """
    counts = [0] * len(regions)
    for idx, row in enumerate(regions):
        for entry in vcf.fetch(row[0], row[1], row[2]):
            if within:
                ent_start, ent_end = truvari.entry_boundaries(entry)
                if not (row[1] <= ent_start and ent_end <= row[2]):
                    continue
            counts[idx] += 1
    return counts

def benchdir_count_entries(benchdir, regions, within=False):
    """
    Count the number of variants per bed region in Truvari bench directory by state

    Returns a pandas dataframe of the counts
    """
    vcfs = {'tpbase': pysam.VariantFile(os.path.join(benchdir, 'tp-base.vcf.gz')),
            'tp': pysam.VariantFile(os.path.join(benchdir, 'tp-comp.vcf.gz')),
            'fn': pysam.VariantFile(os.path.join(benchdir, 'fn.vcf.gz')),
            'fp': pysam.VariantFile(os.path.join(benchdir, 'fp.vcf.gz'))
            }
    data = {}
    for name, vcf in vcfs.items():
        data[name] = count_entries(vcf, regions, within)
    return pd.DataFrame(data)

def stratify_main(cmdargs):
    """
    stratify
    """
    args = parse_args(cmdargs)
    read_header = 0 if args.header else None
    regions = pd.read_csv(args.regions, sep='\t', header=read_header)
    r_list = regions.to_numpy().tolist() # the methods expect lists
    if os.path.isdir(args.in_vcf):
        counts = benchdir_count_entries(args.in_vcf, r_list, args.within)[["tpbase", "tp", "fn", "fp"]]
    else:
        counts = count_entries(pysam.VariantFile(args.in_vcf), r_list, args.within)
        counts = pd.Series(counts, name="count").to_frame()
    counts.index = regions.index
    regions = regions.join(counts)
    regions.to_csv(args.output, header=args.header, index=False, sep='\t')
