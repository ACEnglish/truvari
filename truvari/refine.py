"""
Automated Truvari bench result refinement with phab
"""
import os
import sys
import json
import shutil
import logging
import argparse
from argparse import Namespace

import pandas as pd

from pysam import bcftools
from intervaltree import IntervalTree

import truvari

def read_json(fn):
    """
    Parse json and return dict
    """
    ret = None
    with open(fn, 'r') as fh:
        ret = json.load(fh)
    return ret

def intersect_beds(bed_a, bed_b):
    """
    Return bed_a regions that intersect bed_b regions
    """
    shared = {}
    count = 0
    for chrom in bed_a:
        s_inc = []
        for i in bed_a[chrom]:
            if bed_b[chrom].overlaps(i):
                count += 1
                s_inc.append(i)
        shared[chrom] = IntervalTree(s_inc)
    return shared, count

def tree_to_regions(trees):
    """
    turn an anno tree into a list of lists of chrom start end
    """
    ret = []
    for chrom in trees:
        for intv in trees[chrom]:
            ret.append([chrom, intv.begin, intv.end])
    return ret

def resolve_regions(params, args):
    """
    Figures out or creates the regions we'll be analyzing
    """
    if params["includebed"] is None and args.regions is None:
        logging.error("Bench output didn't use `--includebed` and `--regions` not provided")
        logging.error("Unable to run refine")
        sys.exit(1)
    elif args.regions is None:
        reeval_trees, new_count = truvari.build_anno_tree(params["includebed"], idxfmt="")
        logging.info("Evaluating %d regions", new_count)
    elif args.regions is not None and params["includebed"] is not None:
        a_trees, regi_count = truvari.build_anno_tree(args.regions, idxfmt="")
        b_trees, orig_count = truvari.build_anno_tree(params["includebed"], idxfmt="")
        if args.use_includebed:
            reeval_trees, new_count = intersect_beds(b_trees, a_trees)
            logging.info("%d --includebed reduced to %d after intersecting with %d from --regions",
                     orig_count, new_count, regi_count)

        else:
            reeval_trees, new_count = intersect_beds(a_trees, b_trees)
            logging.info("%d --regions reduced to %d after intersecting with %d from --includebed",
                     regi_count, new_count, orig_count)
    else:
        reeval_trees, count = truvari.build_anno_tree(args.regions, idxfmt="")
        logging.info("%d --regions loaded", count)
    # might need to merge overlaps.
    return reeval_trees


def consolidate_bench_vcfs(benchdir):
    """
    Pull and consolidate base/comp variants from their regions
    """
    bout_name = truvari.make_temp_filename(suffix=".vcf")
    output = bcftools.concat("--no-version", *[os.path.join(benchdir, "tp-base.vcf.gz"),
                            os.path.join(benchdir, "fn.vcf.gz")])
    with open(bout_name, 'w') as fout:
        fout.write(output)
    truvari.compress_index_vcf(bout_name)

    cout_name = truvari.make_temp_filename(suffix=".vcf")
    output = bcftools.concat("--no-version", *[os.path.join(benchdir, "tp-call.vcf.gz"),
                            os.path.join(benchdir, "fp.vcf.gz")])
    with open(cout_name, 'w') as fout:
        fout.write(output)
    truvari.compress_index_vcf(cout_name)
    return bout_name + '.gz', cout_name + '.gz'


def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="refine", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("benchdir", metavar="DIR",
                        help="Truvari bench directory")
    parser.add_argument("-f", "--reference", type=str,
                        help="Indexed fasta used to call variants")
    parser.add_argument("-r", "--regions", default=None,
                        help="Regions to process")
    parser.add_argument("-I", "--use-includebed", action="store_true",
                        help="When intersecting includebed with regions, use includebed coordinates")
    parser.add_argument("-u", "--use-original", action="store_true",
                        help="Use original input VCFs instead of filtered tp/fn/fps")
    parser.add_argument("-t", "--threads", default=1, type=int,
                        help="Number of threads to use (%(default)s)")
    parser.add_argument("-k", "--keep-phab", action="store_true",
                        help="Keep the phab intermediate files")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug, show_version=True)
    return args

def check_params(args):
    """
    Sets up all the parameters from the bench/params.json

    Returns params dict
    """
    param_path = os.path.join(args.benchdir, "params.json")
    if not os.path.exists(param_path):
        logging.error("Bench directory %s doesn't have params.json", param_path)
        sys.exit(1)
    params = read_json(param_path)

    if params["reference"] is None and args.reference is None:
        logging.error("Reference not in params.json or given as a parameter to refine")
        sys.exit(1)
    elif args.reference is None:
        args.reference = params["reference"]

    if not os.path.exists(args.reference):
        logging.error("Reference %s does not exist", args.reference)
        sys.exit(1)

    # Setup prefix
    params["cSample"] = "p:" + params["cSample"]

    phdir = os.path.join(args.benchdir, 'phab')
    if os.path.exists(phdir):
        logging.error("Directory %s exists. Cannot run refine", phdir)
        sys.exit(1)

    bhdir = os.path.join(args.benchdir, 'phab_bench')
    if os.path.exists(bhdir):
        logging.error("Directory %s exists. Cannot run refine", bhdir)
        sys.exit(1)

    # Should check that bench dir has compressed/indexed vcfs for fetching
    check_fail = False
    for i in ["tp-base.vcf.gz", "tp-call.vcf.gz", "fn.vcf.gz", "fp.vcf.gz"]:
        if not os.path.exists(os.path.join(args.benchdir, i)):
            logging.error("Benchdir doesn't have compressed/indexed %s", i)
            check_fail = True
    if check_fail:
        sys.exit(1)

    os.makedirs(phdir)
    return params

def refine_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)

    params = check_params(args)

    reeval_trees = resolve_regions(params, args)
    regions = tree_to_regions(reeval_trees)

    # Stratify.
    counts = truvari.benchdir_count_entries(args.benchdir, regions)[["tpbase", "tp", "fn", "fp"]]
    regions = pd.DataFrame(regions, columns=["chrom", "start", "end"])
    counts.index = regions.index
    counts.columns = ["in_tpbase", "in_tp", "in_fn", "in_fp"]
    regions = regions.join(counts)

    # Figure out which to reevaluate
    regions["refined"] = (regions["in_fn"] > 0) & (regions["in_fp"] > 0)
    logging.info("%d regions to be refined", regions["refined"].sum())

    reeval_bed = truvari.make_temp_filename(suffix=".bed")
    regions[regions["refined"]].to_csv(reeval_bed, sep='\t', header=False, index=False)

    # Set the input VCFs for phab
    if args.use_original:
        base_vcf, comp_vcf = params["base"], params["comp"]
    else:
        base_vcf, comp_vcf = consolidate_bench_vcfs(args.benchdir)

    # Send the vcfs to phab
    to_eval_coords = regions[regions["refined"]][["chrom", "start", "end"]].to_numpy().tolist()
    phab_dir = os.path.join(args.benchdir, "phab")
    truvari.phab_multi(base_vcf, args.reference, phab_dir, to_eval_coords,
                       comp_vcf=comp_vcf, prefix_comp=True, threads=args.threads)

    phab_vcf = os.path.join(args.benchdir, "phab.output.vcf.gz")
    truvari.consolidate_phab_vcfs(phab_dir, phab_vcf)

    # Now run bench on the phab harmonized variants
    m_args = Namespace(**params)
    m_args.no_ref = 'a'
    m_args.output = os.path.join(args.benchdir, "phab_bench")
    m_args.base = phab_vcf
    m_args.comp = phab_vcf
    m_args.includebed = reeval_bed
    truvari.run_bench(m_args)

    # update the output variant counts
    counts = truvari.benchdir_count_entries(m_args.output, to_eval_coords)[["tpbase", "tp", "fn", "fp"]]
    counts.index = regions[regions['refined']].index
    counts.columns = ["out_tpbase", "out_tp", "out_fn", "out_fp"]
    regions = regions.join(counts)
    for i in ["tpbase", "tp", "fn", "fp"]:
        regions[f"out_{i}"] = regions[f"in_{i}"].where(~regions['refined'], regions[f"out_{i}"])

    summary = truvari.StatsBox()
    summary["TP-base"] = int(regions["out_tpbase"].sum())
    summary["TP-call"] = int(regions["out_tp"].sum())
    summary["FP"] = int(regions["out_fp"].sum())
    summary["FN"] = int(regions["out_fn"].sum())
    summary["base cnt"] = summary["TP-base"] + summary["FN"]
    summary["call cnt"] = summary["TP-call"] + summary["FP"]
    # Still don't have genotype checks
    summary.calc_performance()

    with open(os.path.join(args.benchdir, 'refine.summary.json'), 'w') as fout:
        json.dump(summary, fout, indent=4)
    logging.info(json.dumps(summary, indent=4))

    regions.to_csv(os.path.join(args.benchdir, 'refine.counts.txt'), sep='\t', index=False)

    if not args.keep_phab:
        shutil.rmtree(phab_dir)
    logging.info("Finished refine")
