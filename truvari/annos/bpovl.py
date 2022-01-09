"""
Creates intersection of features in an annotation file with SVs' breakpoints and overlap
"""
import gzip
import logging
import argparse
from collections import defaultdict

import pysam
import joblib
import pandas as pd
from intervaltree import IntervalTree

import truvari
# config of [chrom_idx, begin_idx, end_idx, one-based, comment]
PRESET_FMTS = {'bed': [0, 1, 2, False, '#'],
               'gff': [0, 3, 4, True, '#']}


def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="bpovl", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, default="/dev/stdin",
                        help="VCF to annotate (stdin)")
    parser.add_argument("-a", "--anno", type=str, required=True,
                        help="Tab-delimited annotation file")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output joblib DataFrame")
    parser.add_argument("--spanmin", type=int, default=50,
                        help="Minimum span of SVs to annotate (%(default)s)")
    parser.add_argument("--spanmax", type=int, default=50000,
                        help="Maximum span of SVs to annotate (%(default)s)")
    annosg = parser.add_argument_group("Annotation File Arguments")
    annosg.add_argument("-p", "--preset", choices=PRESET_FMTS.keys(), default=None,
                        help=("Annotation format. This option overwrites "
                              "-s, -b, -e, -c and -1 (%(default)s)"))
    annosg.add_argument("-c", "--comment", type=str, default="#",
                        help="Skip lines started with character. (%(default)s)")
    annosg.add_argument("-s", "--sequence", type=int, default=0,
                        help="Column of sequence/chromosome name. (%(default)s)")
    annosg.add_argument("-b", "--begin", type=int, default=1,
                        help="Column of start chromosomal position. (%(default)s)")
    annosg.add_argument("-e", "--end", type=int, default=2,
                        help="Column of end chromosomal position. (%(default)s)")
    # The end column can be the same as the start column. [2]
    annosg.add_argument("-1", "--one-based", action='store_true',
                        help=("The position in the anno file is 1-based "
                              "rather than 0-based. (%(default)s)"))

    args = parser.parse_args(args)
    if args.preset is not None:
        args.anno_psets = PRESET_FMTS[args.preset]
    else:
        args.anno_psets = [args.sequence, args.begin, args.end,
                           args.one_based, args.comment]
    truvari.setup_logging()
    return args


def build_anno_tree(filename, chrom_col=0, start_col=1, end_col=2, one_based=False, comment='#'):
    """
    Build an dictionary of IntervalTrees for each chromosome from tab-delimited annotation file
    """
    def gz_hdlr(fn):
        with gzip.open(fn) as fh:
            for line in fh:
                yield line.decode()

    def fh_hdlr(fn):
        with open(fn) as fh:
            for line in fh:
                yield line

    correction = 1 if one_based else 0
    tree = defaultdict(IntervalTree)
    if filename.endswith('.gz'):
        fh = gz_hdlr(filename)
    else:
        fh = fh_hdlr(filename)

    idx = 0
    for line in fh:
        if line.startswith(comment):
            continue
        data = line.strip().split('\t')
        chrom = data[chrom_col]
        start = int(data[start_col]) - correction
        end = int(data[end_col])
        tree[chrom].addi(start, end, data=idx)
        idx += 1
    return tree, idx


def bpovl_main(cmdargs):
    """
    Main method
    """
    args = parse_args(cmdargs)
    in_vcf = pysam.VariantFile(args.input)
    anno_tree, anno_cnt = build_anno_tree(args.anno, *args.anno_psets)
    out_rows = []
    logging.info("Loaded %d annotations", anno_cnt)
    hit_cnt = 0
    for entry in in_vcf:
        has_hit = False
        start, end = truvari.entry_boundaries(entry)
        span = abs(end - start)
        if span < args.spanmin or span > args.spanmax:
            continue
        key = truvari.entry_to_key(entry)
        for anno_idx in anno_tree[entry.chrom].at(start):
            has_hit = True
            out_rows.append([key, 'start_bnd', anno_idx.data])

        for anno_idx in anno_tree[entry.chrom].at(end):
            has_hit = True
            out_rows.append([key, 'end_bnd', anno_idx.data])

        for anno_idx in anno_tree[entry.chrom].overlap(start, end):
            has_hit = True
            if anno_idx.begin >= start and anno_idx.end <= end:
                out_rows.append([key, 'overlaps', anno_idx.data])
            elif anno_idx.begin < start and anno_idx.end > end:
                out_rows.append([key, 'contains', anno_idx.data])
        hit_cnt += has_hit
    logging.info("%d SVs hit annotations", hit_cnt)
    out = pd.DataFrame(out_rows, columns=["vcf_key",
                                          "intersection",
                                          "anno_key"])
    joblib.dump(out, args.output)
    logging.info("Finished bpovl")
