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
    parser.add_argument("in_vcf", metavar="VCF",
                        help="Compressed Index VCF")
    parser.add_argument("regions", metavar="BED",
                        help="Regions to process")
    parser.add_argument("-o", "--output", metavar="OUT", default="/dev/stdout",
                        help="Output bed-like file")
    parser.add_argument("-w", "--within", action="store_true",
                        help="Only count variants contained completely within region boundaries")
    parser.add_argument("-b", "--bench-dir", action="store_true",
                        help="in_vcf is a bench output directory. Parse each VCF")
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
    counts = []
    for row in regions:
        m_cnt = 0
        for entry in vcf.fetch(row[0], row[1], row[2]):
            if within:
                ent_start, ent_end = truvari.entry_boundaries(entry)
                if not row[1] <= ent_start and ent_end <= row[2]:
                    continue
            m_cnt += 1
        counts.append(m_cnt)
    return counts

def benchdir_count_entries(benchdir, regions, within=False):
    """
    Count the number of variants per bed region in Truvari bench directory by state
    Returns a pandas dataframe of the counts
    """
    vcfs = {'tpbase': pysam.VariantFile(os.path.join(benchdir, 'tp-base.vcf.gz')),
            'tp': pysam.VariantFile(os.path.join(benchdir, 'tp-call.vcf.gz')),
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

    regions = pd.read_csv(args.regions, sep='\t', header=None)
    r_list = regions.to_numpy().tolist() # the methods expect lists
    if args.bench_dir:
        counts = benchdir_count_entries(args.in_vcf, r_list, args.within)[["tpbase", "tp", "fn", "fp"]]
    else:
        counts = count_entries(pysam.VariantFile(args.in_vcf), r_list, args.within)
        counts = pd.Series(counts, name="count").to_frame()
    counts.index = regions.index
    regions = regions.join(counts)
    regions.to_csv(args.output, header=False, index=False, sep='\t')
