"""
Automated Truvari bench result refinement
"""
import os
import sys
import json
import logging
import argparse
from typing import List
from dataclasses import dataclass, field

import pysam
from intervaltree import IntervalTree, Interval

import truvari
from truvari.bench import StatsBox

def read_json(fn):
    """
    Parse json and return dict
    """
    ret = None
    with open(fn, 'r') as fh:
        ret = json.load(fh)
    return ret

def intersect_beds(includebed, regions):
    """
    Remove includebed regions that do not intersect subset regions
    Return the unique and shared regions as dict of IntervalTrees
    """
    unique_include = {}
    shared_include = {}
    for chrom in includebed:
        u_inc = []
        s_inc = []
        for i in includebed[chrom]:
            if not regions[chrom].overlaps(i):
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
@dataclass
class ReevalRegion:
    """
    Class for keeping track of a region for processing.
    """
    name: str
    benchdir: str
    params: dict
    chrom: str
    start: int
    end: int
    needs_reeval: bool = True
    # These may need to become a temporary file's name
    base_calls: List = field(default_factory=lambda: [])
    comp_calls: List = field(default_factory=lambda: [])
    in_tp_base_count: int = 0
    in_tp_call_count: int = 0
    in_fp_count: int = 0
    in_fn_count: int = 0
    out_tp_base_count: int = 0
    out_tp_comp_count: int = 0
    out_fp_count: int = 0
    out_fn_count: int = 0
    
    def update_summary(self, summary):
        """
        Given a StatsBox (dict), update the counts based on what we've observed

        This is currently done stupidly. I'm sure there's a cleaner logic for updating the counts, 
        but I'll refactor it later
        """
        summary["TP-base"] -= self.in_tp_base_count
        summary["TP-call"] -= self.in_tp_call_count
        summary["FN"] -= self.in_fn_count
        summary["FP"] -= self.in_fp_count
        summary["base cnt"] -= self.in_tp_base_count + self.in_fn_count
        summary["call cnt"] -= self.in_tp_call_count + self.in_fp_count

        summary["TP-base"] += self.out_tp_base_count
        summary["TP-call"] += self.out_tp_comp_count
        summary["FN"] += self.out_fn_count
        summary["FP"] += self.out_fp_count
        summary["base cnt"] += self.out_tp_base_count + self.out_fn_count
        summary["call cnt"] += self.out_tp_comp_count + self.out_fp_count

def region_needs_reeval(region):
    """
    Determines if a ReevalRegion even needs to be re-evaluated

    Note: eventually I would like this method to also determine which EVALS method is best.
    When we implement that, we'll add a new member to the data class and... I might can
    design this in now...
    """
    region.needs_reeval = region.in_fn_count > 0 and region.in_fp_count > 0
    logging.debug("region %s needs reeval? %s", region.name, region.needs_reeval)
    return region

def fetch_originals(region):
    """
    Grabs the variants that are needed for processing

    Populates the RevalRegion's `base_calls` and `comp_calls`
    """
    fetch_bench_vcfs(region, populate_calls=False)
    region.base_calls = [_ for _ in pysam.VariantFile(region.params["base"]).fetch(region.chrom, region.start, region.end)]
    region.comp_calls = [_ for _ in pysam.VariantFile(region.params["comp"]).fetch(region.chrom, region.start, region.end)]
    return region   

def fetch_bench_vcfs(region, populate_calls=True):
    """
    Grabs the variants that are needed for processing
    Populate the input state counts from truvari vcfs

    Populates the RevalRegion's `base_calls` and `comp_calls`
    """
    vcf = pysam.VariantFile(os.path.join(region.benchdir, "tp-base.vcf.gz"))
    for entry in vcf.fetch(region.chrom, region.start, region.end):
        region.in_tp_base_count += 1
        if populate_calls:
            region.base_calls.append(entry)

    vcf = pysam.VariantFile(os.path.join(region.benchdir, "fn.vcf.gz"))
    for entry in vcf.fetch(region.chrom, region.start, region.end):
        region.in_fn_count += 1
        if populate_calls:
            region.base_calls.append(entry)

    vcf = pysam.VariantFile(os.path.join(region.benchdir, "tp-call.vcf.gz"))
    for entry in vcf.fetch(region.chrom, region.start, region.end):
        region.in_tp_call_count += 1
        if populate_calls:
            region.comp_calls.append(entry)

    vcf = pysam.VariantFile(os.path.join(region.benchdir, "fp.vcf.gz"))
    for entry in vcf.fetch(region.chrom, region.start, region.end):
        region.in_fp_count += 1
        if populate_calls:
            region.comp_calls.append(entry)

    return region

def phab_eval(region):
    """
    Analyze a ReevalRegion with phab
    """
    # Short circuit
    if not region.needs_reeval:
        return region
    # Just move one over for fun
    region.out_fn_count = region.in_fn_count - 1
    region.out_fp_count = region.in_fp_count - 1
    region.out_tp_base_count = region.in_tp_base_count + 1
    region.out_tp_comp_count = region.in_tp_comp-count + 1

    return region
    #raise NotImplementedError("In progress")

def hap_eval(region):
    """
    Analyze a ReevalRegion with hap_eval
    """
    raise NotImplementedError("In progress")
    # Short circuit
    if not region.needs_reeval:
        return region

EVALS = {'phab':phab_eval, 'hap':hap_eval}

def rebench(benchdir, params, summary, reeval_trees, use_original=False, eval_method='phab', workers=1):
    """
    Sets up inputs before calling specified EVAL method
    Will take the output from the EVAL method to then update summary

    :param `benchdir`: Truvari bench output directory
    :type `benchdir`: :class:`str`
    :param `params`: bench params.json
    :type `benchdir`: :class:`dict`
    :param `summary`: bench summary.json
    :type `bench`: :class:`dict`
    :param `reeval_trees`: dict of IntervalTrees of regions to be re-evaluated
    :type `reeval_trees`: `dict`
    :param `use_original`: true if use original VCFs, else parse calls in benchdir VCFs
    :type `use_original`: `bool`
    :param `eval_method`:
    :param `workers`: number of workers to use
    :type `workers`: `int`
    """
    # This will process multiple regions, so think about how to separate region parsing from eval
    # If I keep it clean we'll have an opportunity to parallelize. Can maybe make use of callbacks, also
    
    # Build a Region datastructure. These are what we use to pass inputs/output
    data = []
    for chrom in reeval_trees:
        for intv in reeval_trees[chrom]:
            name = f"{chrom}:{intv.begin}-{intv.end}"
            data.append(ReevalRegion(name, benchdir, params, chrom, intv.begin, intv.end))
    
    # Build the pipeline
    pipeline = []

    # Setup variants
    if use_original:
        pipeline.append((fetch_originals))
    else:
        pipeline.append((fetch_bench_vcfs))

    # Is this even a region worth considering?
    pipeline.append((region_needs_reeval))

    # Call the eval_method
    pipeline.append((EVALS[eval_method]))

    # Collect eval_method's results
    for result in truvari.fchain(pipeline, data, workers=workers):
        result.update_summary(summary)

    # Doesn't need to return anything since summary is updated in-place
    return

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="rebench", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("benchdir", metavar="DIR",
                        help="Truvari bench directory")
    parser.add_argument("-f", "--reference", type=str, required=True,
                        help="Indexed fasta used to call variants")
    parser.add_argument("-r", "--regions", default=None,
                        help="Regions to process")
    parser.add_argument("-e", "--eval", default='phab', choices=EVALS.keys(),
                        help="Evaluation procedure (%(default)s)")
    parser.add_argument("-u", "--use-original", action="store_true",
                        help="Use original input VCFs instead of filtered tp/fn/fps")
    parser.add_argument("-t", "--threads", default=1, type=int,
                        help="Number of threads to use (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    return args

def rebench_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)
    param_path = os.path.join(args.benchdir, "params.json")
    if not os.path.exists(param_path):
        logging.error("Bench directory %s doesn't have params.json", param_path)
        sys.exit(1)
    params = read_json(param_path)

    summary_path = os.path.join(args.benchdir, "summary.json")
    if not os.path.exists(summary_path):
        logging.error("Bench directory %s doesn't have summary.json", param_path)
        sys.exit(1)
    summary = StatsBox()
    summary.update(read_json(summary_path))

    if params["includebed"] is None and args.regions is None:
        logging.error("Bench output didn't use `--includebed` and `--regions` not provided")
        logging.error("Unable to run rebench")
        sys.exit(1)

    if params["includebed"] is None:
        btree, _ = truvari.build_anno_tree(args.regions, idxfmt="")
    else:
        btree, _ = truvari.build_anno_tree(params["includebed"], idxfmt="")

    # Logic gap here. If there is no includebed, we don't build fn_trees
    # I think that's okay?...
    if args.regions:
        stree, _ = truvari.build_anno_tree(args.regions)
        fn_trees, reeval_trees = intersect_beds(btree, stree)
    else:
        fn_trees = {}
        reeval_trees = btree

    # Should check that bench dir has compressed/indexed vcfs for fetching
    update_fns(args.benchdir, fn_trees, summary)

    # Will eventually need to pass args for phab and|or hap-eval
    rebench(args.benchdir, params, summary, reeval_trees, args.use_original, args.eval, args.threads)
    # Are we going to annotate the regions with the tp/fp counts and changes.. maybe eventually
    summary.calc_performance()
    print(json.dumps(summary, indent=4))
