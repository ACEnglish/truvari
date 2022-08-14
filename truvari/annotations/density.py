"""
Identify 'dense' and 'sparse' variant windows of the genome
"""
import logging
import argparse
from collections import Counter, defaultdict
import pysam
import joblib
import pandas as pd
from intervaltree import IntervalTree

import truvari

def parse_args(args):
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(prog="density", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-g", "--genome", type=str, required=True,
                        help="Genome bed file")
    parser.add_argument("-i", "--input", type=str, default="/dev/stdin",
                        help="Input VCF (%(default)s)")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output joblib DataFrame")
    parser.add_argument("-m", "--mask", type=str,
                        help="Mask bed file")
    parser.add_argument("-w", "--windowsize", type=truvari.restricted_int, default=10000,
                        help="Window size (%(default)s)")
    parser.add_argument("-t", "--threshold", type=float, default=3,
                        help="std for identifying 'dense' regions (%(default)s)")
    args = parser.parse_args(args)
    truvari.setup_logging()
    return args

def density_main(args):
    """
    Main
    """
    args = parse_args(args)
    tree, cnt = truvari.build_anno_tree(args.genome)
    logging.info("Loaded %d seqs from genome", cnt)

    # mask
    if args.mask:
        mask_cnt = 0
        with open(args.mask) as fh:
            for line in fh:
                mask_cnt += 1
                data = line.strip().split('\t')
                chrom = data[0]
                start = int(data[1])
                end = int(data[2])
                tree[chrom].chop(start, end)
        logging.info("Masked %d regions", mask_cnt)

    # setting new indexes after masking
    new_tree = defaultdict(IntervalTree)
    cnt = 0
    for chrom in tree:
        for intv in tree[chrom]:
            for i in range(intv.begin, intv.end, args.windowsize):
                new_tree[chrom].addi(i, min(intv.end, i + args.windowsize), data=cnt)
                cnt += 1
    logging.info("Made %d %dbp windows", cnt, args.windowsize)
    tree = new_tree

    # Counting
    counts = Counter()
    v = pysam.VariantFile(args.input)
    cnt = 0
    for entry in v:
        for intv in tree[entry.chrom].overlap(entry.start, entry.stop):
            cnt += 1
            counts[intv.data] += 1
    logging.info("Intersected %d variants", cnt)

    # Output
    data = []
    for chrom in tree:
        for intv in tree[chrom]:
            data.append([chrom, intv.begin, intv.end, counts[intv.data]])
    data = pd.DataFrame(data, columns=['chrom', 'start', 'end', 'count'])
    # Do the hotspot work
    desc = data["count"].describe()
    logging.info("Summary\n%s", str(desc))
    hs_threshold = desc["mean"] + (args.threshold * desc["std"])
    logging.info("Setting threshold %f * SD = %f", args.threshold, hs_threshold)

    data["anno"] = None
    data.loc[data["count"] == 0, "anno"] = "sparse"
    data.loc[data["count"] > hs_threshold, "anno"] = "dense"
    logging.info("Density Counts\n%s", str(data["anno"].value_counts(dropna=False)))

    joblib.dump(data, args.output)
