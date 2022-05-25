"""
Creates intersection of features in an annotation file with SVs' breakpoints and overlap
"""
import logging
import argparse

import pysam
import joblib
import pandas as pd

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
    parser.add_argument("--sizemin", type=truvari.restricted_int, default=50,
                        help="Minimum size of variant to annotate (%(default)s)")
    parser.add_argument("--spanmax", type=truvari.restricted_int, default=50000,
                        help="Maximum span of SVs to annotate (%(default)s)")
    annosg = parser.add_argument_group("Annotation File Arguments")
    annosg.add_argument("-p", "--preset", choices=PRESET_FMTS.keys(), default=None,
                        help=("Annotation format. This option overwrites "
                              "-s, -b, -e, -c and -1 (%(default)s)"))
    annosg.add_argument("-c", "--comment", type=str, default="#",
                        help="Skip lines started with character. (%(default)s)")
    annosg.add_argument("-s", "--sequence", type=truvari.restricted_int, default=0,
                        help="Column of sequence/chromosome name. (%(default)s)")
    annosg.add_argument("-b", "--begin", type=truvari.restricted_int, default=1,
                        help="Column of start chromosomal position. (%(default)s)")
    annosg.add_argument("-e", "--end", type=truvari.restricted_int, default=2,
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


def bpovl_main(cmdargs):
    """
    Main method
    """
    args = parse_args(cmdargs)
    in_vcf = pysam.VariantFile(args.input)
    anno_tree, anno_cnt = truvari.build_anno_tree(args.anno, *args.anno_psets)
    out_rows = []
    logging.info("Loaded %d annotations", anno_cnt)
    hit_cnt = 0
    for entry in in_vcf:
        has_hit = False

        start, end = truvari.entry_boundaries(entry)
        span = abs(end - start)
        if span > args.spanmax:
            continue
        svlen = truvari.entry_size(entry)
        if svlen < args.sizemin:
            continue

        key = truvari.entry_to_key(entry)
        for anno_idx in anno_tree[entry.chrom].at(start):
            has_hit = True
            out_rows.append([key, 'start_bnd', anno_idx.data])

        for anno_idx in anno_tree[entry.chrom].at(end):
            has_hit = True
            out_rows.append([key, 'end_bnd', anno_idx.data])

        for anno_idx in anno_tree[entry.chrom].overlap(start, end):
            if start <= anno_idx.begin and anno_idx.end <= end:
                has_hit = True
                out_rows.append([key, 'overlaps', anno_idx.data])
            elif anno_idx.begin <= start and end <= anno_idx.end:
                has_hit = True
                out_rows.append([key, 'contains', anno_idx.data])
        hit_cnt += has_hit
    logging.info("%d SVs hit annotations", hit_cnt)
    out = pd.DataFrame(out_rows, columns=["vcf_key",
                                          "intersection",
                                          "anno_key"])
    joblib.dump(out, args.output)
    logging.info("Finished bpovl")
