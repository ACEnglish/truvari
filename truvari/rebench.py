"""
Automated Truvari bench result refinement
"""
import os
import json
import argparse

import pysam
from intervaltree import IntervalTree

import truvari
from truvari.bench import StatsBox

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="rebench", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("benchdir", metavar="DIR",
                        help="Truvari bench directory")
    parser.add_argument("--subset", default=None,
                        help="Subset of regions to process")
    args = parser.parse_args(args)
    return args

def read_json(fn):
    """
    Parse and return a json
    """
    ret = None
    with open(fn, 'r') as fh:
        ret = json.load(fh)
    return ret


def intersect_beds(includebed, subset):
    """
    Remove includebed regions that do not intersect subset regions
    Return the unique and shared bed files (as dict of IntervalTrees)
    """
    unique_include = {}
    shared_include = {}
    for chrom in includebed:
        u_inc = []
        s_inc = []
        for i in includebed[chrom]:
            if not subset[chrom].overlaps(i):
                u_inc.append(i)
            else:
                s_inc.append(i)
        unique_include[chrom] = IntervalTree(u_inc)
        shared_include[chrom] = IntervalTree(s_inc)
    return unique_include, shared_include

def update_fns(benchdir, fn_trees, summary):
    """
    For all FP calls in the fn_trees regions, update the FN counts
    summary is updated in-place
    """
    vcf = pysam.VariantFile(os.path.join(benchdir, "fn.vcf.gz"))
    for chrom in fn_trees:
        for intv in fn_trees[chrom]:
            for _ in vcf.fetch(chrom, intv.begin, intv.end):
                summary["FN"] -= 1
                summary["base cnt"] -= 1
                # I think I can use regionvcfiter to make the output rebench.fn.vcf  more easily

#def recompare(*args, **kwargs):
# """
# 1. Need to hook in phab first.
#   truvari.phab
#   There's a problem here
# 2. Then we need a way to parse each of the phab results
# 3. We can have a flow control here that does hap-eval instead
# """

def rebench_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)
    params = read_json(os.path.join(args.benchdir, "params.json"))
    summary = StatsBox()
    summary.update(read_json(os.path.join(args.benchdir, "summary.json")))

    print(params["includebed"])
    # Should check that its used. If not, this is an error
    btree, _ = truvari.build_anno_tree(params["includebed"], idxfmt="")
    if args.subset:
        stree, _ = truvari.build_anno_tree(args.subset)
        fn_trees, reeval_trees = intersect_beds(btree, stree)
    else:
        fn_trees = {}
        reeval_trees = btree

    # Should check that bench dir has compressed/indexed vcfs for fetching
    update_fns(args.benchdir, fn_trees, summary)

    # Will eventually need to pass args for phab and|or hap-eval
    #recompare(args.benchdir, reeval_trees, summary)
    print(type(reeval_trees))
    summary.calc_performance()
    print(json.dumps(summary, indent=4))
