"""
Count variants per-region in vcf
"""
import os
import logging
import argparse
import multiprocessing
from functools import partial
from collections import defaultdict

import pysam
import numpy as np
import pandas as pd
from intervaltree import IntervalTree

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


def count_entries(vcf, chroms, regions, within):
    """
    Count the number of variants per bed region a VCF
    regions is a list of lists with sub-lists having 0:3 of chrom, start, end
    Returns a list of the counts in the same order as the regions
    """
    if isinstance(vcf, str):
        vcf = pysam.VariantFile(vcf)
    tree = defaultdict(IntervalTree)
    counts_idx = {}
    counts = [0] * len(regions)
    for idx, row in enumerate(zip(chroms, regions)):
        chrom, coords = row
        start, end = coords
        end += 1
        tree[chrom].addi(start, end)
        counts_idx[(chrom, start, end)] = idx
    for _, location in truvari.region_filter(vcf, tree, within, True):
        key = (location[0], location[1].begin, location[1].end)
        counts[counts_idx[key]] += 1
    return counts


def benchdir_count_entries(benchdir, regions, within=False, threads=4):
    """
    Count the number of variants per bed region in Truvari bench directory by state

    Returns a pandas dataframe of the counts
    """
    names = ['tpbase', 'tp', 'fn', 'fp']
    vcfs = [os.path.join(benchdir, 'tp-base.vcf.gz'),
            os.path.join(benchdir, 'tp-comp.vcf.gz'),
            os.path.join(benchdir, 'fn.vcf.gz'),
            os.path.join(benchdir, 'fp.vcf.gz')]
    if isinstance(regions, pd.DataFrame):
        regions = regions.to_numpy().tolist()  # the methods expect lists
    chroms = np.array([_[0] for _ in regions])
    intvs = np.array([[_[1], _[2]] for _ in regions])
    method = partial(count_entries, chroms=chroms,
                     regions=intvs, within=within)
    data = {}
    with multiprocessing.Pool(threads) as pool:
        for name, counts in zip(names, pool.map(method, vcfs)):
            data[name] = counts
    return pd.DataFrame(data)


def stratify_main(cmdargs):
    """
    stratify
    """
    args = parse_args(cmdargs)
    read_header = 0 if args.header else None
    regions = pd.read_csv(args.regions, sep='\t', header=read_header)
    r_list = regions.to_numpy().tolist()  # the methods expect lists
    if os.path.isdir(args.in_vcf):
        counts = benchdir_count_entries(args.in_vcf, r_list, args.within)
    else:
        chroms = np.array([_[0] for _ in r_list])
        intvs = np.array([[_[1], _[2]] for _ in r_list])
        counts = count_entries(pysam.VariantFile(
            args.in_vcf), chroms, intvs, args.within)
        counts = pd.Series(counts, name="count").to_frame()
    counts.index = regions.index
    regions = regions.join(counts)
    regions.to_csv(args.output, header=args.header, index=False, sep='\t')
    logging.info("Finished")
