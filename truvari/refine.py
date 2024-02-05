"""
Automated Truvari bench result refinement with phab
"""
import os
import sys
import json
import logging
import argparse
from functools import partial
import multiprocessing as mp
from argparse import Namespace
from collections import defaultdict

import pysam
import pandas as pd
from pysam import bcftools
from intervaltree import IntervalTree

import truvari
from truvari.phab import check_requirements as phab_check_requirements
from truvari.phab import DEFAULT_MAFFT_PARAM

PHAB_BUFFER = 100  # hard-coded parameter :(


def intersect_beds(bed_a, bed_b):
    """
    Return bed_a regions that intersect bed_b regions
    """
    shared = {}
    count = 0
    for chrom in bed_a:
        a_intv = sorted(bed_a[chrom])
        bed_b[chrom].merge_overlaps()
        b_intv = sorted(bed_b[chrom])
        a_idx = 0
        b_idx = 0
        s_inc = []
        while a_idx < len(a_intv) and b_idx < len(b_intv):
            if a_intv[a_idx].end < b_intv[b_idx].begin:
                a_idx += 1
                continue
            if b_intv[b_idx].end < a_intv[a_idx].begin:
                b_idx += 1
                continue
            # not before and not after means overlap
            count += 1
            s_inc.append(a_intv[a_idx])
            a_idx += 1
        shared[chrom] = IntervalTree(s_inc)
    return shared, count


def resolve_regions(params, args):
    """
    Figures the regions we'll be analyzing
    """
    logging.info("Setting up regions")
    if params.includebed is None and args.regions is None:
        logging.error(
            "Bench output didn't use `--includebed` and `--regions` not provided")
        logging.error("Unable to run refine")
        sys.exit(1)
    elif args.regions is None:
        reeval_trees, new_count = truvari.build_anno_tree(
            params.includebed, idxfmt="")
        logging.info("Evaluating %d regions", new_count)
    elif args.regions is not None and params.includebed is not None:
        a_trees, regi_count = truvari.build_anno_tree(args.regions, idxfmt="")
        b_trees, orig_count = truvari.build_anno_tree(
            params.includebed, idxfmt="")
        if args.use_region_coords:
            reeval_trees, new_count = intersect_beds(a_trees, b_trees)
            logging.info("%d --regions reduced to %d after intersecting with %d from --includebed",
                         regi_count, new_count, orig_count)
        else:
            reeval_trees, new_count = intersect_beds(b_trees, a_trees)
            logging.info("%d --includebed reduced to %d after intersecting with %d from --regions",
                         orig_count, new_count, regi_count)
    else:
        reeval_trees, count = truvari.build_anno_tree(args.regions, idxfmt="")
        logging.info("%d --regions loaded", count)

    return [[chrom, intv.begin, intv.end - 1]
            for chrom, all_intv in reeval_trees.items()
            for intv in sorted(all_intv)]


def consolidate_bench_vcfs(benchdir):
    """
    Pull and consolidate base/comp variants from their regions
    """
    bout_name = truvari.make_temp_filename(suffix=".vcf")
    output = bcftools.concat("--no-version", "-a", os.path.join(benchdir, "tp-base.vcf.gz"),
                             os.path.join(benchdir, "fn.vcf.gz"))
    with open(bout_name, 'w') as fout:
        fout.write(output)
    truvari.compress_index_vcf(bout_name)

    cout_name = truvari.make_temp_filename(suffix=".vcf")
    output = bcftools.concat("--no-version", "-a", os.path.join(benchdir, "tp-comp.vcf.gz"),
                             os.path.join(benchdir, "fp.vcf.gz"))
    with open(cout_name, 'w') as fout:
        fout.write(output)
    truvari.compress_index_vcf(cout_name)
    return bout_name + '.gz', cout_name + '.gz'


def initial_stratify(benchdir, regions, threads=1):
    """
    Get the stratify counts from the initial bench results
    """
    counts = truvari.benchdir_count_entries(benchdir, regions, True, threads)
    regions = pd.DataFrame(regions, columns=["chrom", "start", "end"])
    counts.index = regions.index
    counts.columns = ["in_tpbase", "in_tp", "in_fn", "in_fp"]
    regions = regions.join(counts)
    return regions


def original_stratify(base_vcf, comp_vcf, regions):
    """
    Get the variant counts from the original vcfs and return a list of refine bools
    The logic here is that `(regions["in_fn"] > 0) | (regions["in_fp"] > 0)` may create
    a number of regions which don't benefit from phab because either is without variants.
    So, for that subset we'll double check that base AND comp have any variants to compare.
    Some regions will still be evaluated for no good reason e.g. PASS only or a rogue snp.
    """
    to_eval = (regions["in_fn"] > 0) | (regions["in_fp"] > 0)
    candidates = regions[to_eval]
    chroms = candidates["chrom"].to_numpy()
    intvs = candidates[["start", "end"]].to_numpy()
    method = partial(truvari.count_entries, chroms=chroms,
                     regions=intvs, within=True)
    results = []
    with mp.Pool(2, maxtasksperchild=1) as pool:
        for result in pool.map(method, [base_vcf, comp_vcf]):
            results.append(pd.Series(result, index=candidates.index))
    return to_eval & (results[0] != 0) & (results[1] != 0)


def refined_stratify(outdir, to_eval_coords, regions, threads=1):
    """
    update regions in-place with the output variant counts
    """
    counts = truvari.benchdir_count_entries(
        outdir, to_eval_coords, True, threads)
    counts.index = regions[regions['refined']].index
    counts.columns = ["out_tpbase", "out_tp", "out_fn", "out_fp"]
    regions = regions.join(counts)
    for i in ["tpbase", "tp", "fn", "fp"]:
        regions[f"out_{i}"] = regions[f"in_{i}"].where(
            ~regions['refined'], regions[f"out_{i}"])
    return regions


def make_variant_report(regions):
    """
    Count all the variants, log and write to output
    """
    summary = truvari.StatsBox()
    summary["TP-base"] = int(regions["out_tpbase"].sum())
    summary["TP-comp"] = int(regions["out_tp"].sum())
    summary["FP"] = int(regions["out_fp"].sum())
    summary["FN"] = int(regions["out_fn"].sum())
    summary["base cnt"] = summary["TP-base"] + summary["FN"]
    summary["comp cnt"] = summary["TP-comp"] + summary["FP"]
    summary.calc_performance()
    return summary


def recount_variant_report(orig_dir, phab_dir, regions):
    """
    Count original variants not in refined regions and
    consolidate with the refined counts.
    """
    tree = defaultdict(IntervalTree)
    n_regions = regions[regions["refined"]].copy()
    for _, row in n_regions.iterrows():
        tree[row['chrom']].addi(row['start'], row['end'] + 1)

    summary = truvari.StatsBox()
    with open(os.path.join(phab_dir, "summary.json")) as fh:
        summary.update(json.load(fh))

    # Adding the original counts to the updated phab counts
    vcf = pysam.VariantFile(os.path.join(orig_dir, 'tp-base.vcf.gz'))
    tpb = len(list(truvari.region_filter(vcf, tree, False)))
    summary["TP-base"] += tpb

    vcf = pysam.VariantFile(os.path.join(orig_dir, 'tp-comp.vcf.gz'))
    tpc = len(list(truvari.region_filter(vcf, tree, False)))
    summary["TP-comp"] += tpc

    vcf = pysam.VariantFile(os.path.join(orig_dir, 'fp.vcf.gz'))
    fp = len(list(truvari.region_filter(vcf, tree, False)))
    summary["FP"] += fp

    vcf = pysam.VariantFile(os.path.join(orig_dir, 'fn.vcf.gz'))
    fn = len(list(truvari.region_filter(vcf, tree, False)))
    summary["FN"] += fn

    summary["comp cnt"] += tpc + fp
    summary["base cnt"] += tpb + fn
    summary.calc_performance()
    return summary


def make_region_report(data):
    """
    Given a refine counts DataFrame, calculate the performance of
    PPV, TNR, etc. Also adds 'state' column to regions inplace
    """
    false_pos = data['out_fp'] != 0
    false_neg = data['out_fn'] != 0
    any_false = false_pos | false_neg

    true_positives = (data['out_tp'] != 0) & (
        data['out_tpbase'] != 0) & ~any_false

    true_negatives = (
        data[['out_tpbase', 'out_tp', 'out_fn', 'out_fp']] == 0).all(axis=1)

    state_map = defaultdict(lambda: 'UNK')
    state_map.update({(True, False, False, False): 'TP',
                      (False, True, False, False): 'TN',
                      (False, False, True, False): 'FP',
                      (False, False, False, True): 'FN',
                      (False, False, True, True): 'FN,FP'})

    data['state'] = [state_map[_]
                     for _ in zip(true_positives, true_negatives, false_pos, false_neg)]

    tot = (true_positives + false_pos + false_neg + true_negatives).astype(int)
    und = data.iloc[tot[tot != 1].index]

    result = {}

    base = (data['out_tpbase'] != 0) | (data['out_fn'] != 0)
    baseP = int(base.sum())
    baseN = int((~base).sum())
    comp = (data['out_tp'] != 0) | (data['out_fp'] != 0)
    compP = int(comp.sum())
    compN = int((~comp).sum())

    result["TP"] = int(true_positives.sum())
    result["TN"] = int(true_negatives.sum())
    result["FP"] = int(false_pos.sum())
    result["FN"] = int(false_neg.sum())
    result["base P"] = baseP
    result["base N"] = baseN
    result["comp P"] = compP
    result["comp N"] = compN
    # precision
    result["PPV"] = result["TP"] / \
        result["comp P"] if result["comp P"] != 0 else None
    # recall
    result["TPR"] = result["TP"] / \
        result["base P"] if result["base P"] != 0 else None
    # specificity
    result["TNR"] = result["TN"] / \
        result["base N"] if result["base N"] != 0 else None
    # negative predictive value
    result["NPV"] = result["TN"] / \
        result["comp N"] if result["comp N"] != 0 else None
    # accuracy
    if result["base P"] + result["base N"] != 0:
        result["ACC"] = (result["TP"] + result["TN"]) / \
            (result["base P"] + result["base N"])
    else:
        result["ACC"] = None
    if result["TPR"] is not None and result["TNR"] is not None:
        result["BA"] = (result["TPR"] + result["TNR"]) / 2
    else:
        result["BA"] = None

    if result["PPV"] and result["TPR"]:
        result["F1"] = 2 * ((result["PPV"] * result["TPR"]) /
                            (result["PPV"] + result["TPR"]))
    else:
        result["F1"] = None
    result["UND"] = len(und)

    return result


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
    parser.add_argument("-u", "--use-original-vcfs", action="store_true",
                        help="Use original input VCFs instead of filtered tp/fn/fps")
    parser.add_argument("-U", "--use-region-coords", action="store_true",
                        help="When subsetting the includebed with regions, use region coordinates")
    parser.add_argument("-R", "--recount", action="store_true",
                        help=("Recount all variants for refine.variant_summary instead of only "
                              "those in analyzed regions"))
    parser.add_argument("-t", "--threads", default=4, type=int,
                        help="Number of threads to use (%(default)s)")
    parser.add_argument("-a", "--align", type=str, choices=["mafft", "wfa", "poa"], default="mafft",
                        help="Alignment method for phab (%(default)s)")
    parser.add_argument("-m", "--mafft-params", type=str, default=DEFAULT_MAFFT_PARAM,
                        help="Parameters for mafft, wrap in a single quote (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    return args


def check_params(args):
    """
    Sets up all the parameters from the bench/params.json

    Returns as a Namespace
    """
    check_fail = False
    param_path = os.path.join(args.benchdir, "params.json")
    if not os.path.exists(param_path):
        check_fail = True
        logging.error(
            "Bench directory %s doesn't have params.json", param_path)

    params = None
    with open(param_path, 'r') as fh:
        params = json.load(fh)

    p_file = os.path.join(args.benchdir, "phab.output.vcf.gz")
    if os.path.exists(p_file):
        check_fail = True
        logging.error("Phab output file %s already exists", p_file)

    if params["reference"] is None and args.reference is None:
        check_fail = True
        logging.error(
            "Reference not in params.json or given as a parameter to refine")
    elif args.reference is None:
        args.reference = params["reference"]

    if not os.path.exists(args.reference):
        check_fail = True
        logging.error("Reference %s does not exist", args.reference)

    # Setup prefix
    params["cSample"] = "p:" + params["cSample"]

    bhdir = os.path.join(args.benchdir, 'phab_bench')
    if os.path.exists(bhdir):
        check_fail = True
        logging.error("Directory %s exists. Cannot run refine", bhdir)

    # Should check that bench dir has compressed/indexed vcfs for fetching
    for i in ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fn.vcf.gz", "fp.vcf.gz"]:
        if not os.path.exists(os.path.join(args.benchdir, i)):
            check_fail = True
            logging.error("Benchdir doesn't have compressed/indexed %s", i)

    if check_fail:
        logging.error("Please fix parameters")
        sys.exit(1)

    return Namespace(**params)


def refine_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)
    if phab_check_requirements(args.align):
        logging.error("Couldn't run Truvari. Please fix parameters\n")
        sys.exit(100)

    params = check_params(args)
    truvari.setup_logging(args.debug,
                          truvari.LogFileStderr(os.path.join(
                              args.benchdir, "refine.log.txt")),
                          show_version=True)
    logging.info("Params:\n%s", json.dumps(vars(args), indent=4))

    # Stratify.
    regions = initial_stratify(args.benchdir,
                               resolve_regions(params, args),
                               args.threads)

    # Figure out which to reevaluate
    if args.use_original_vcfs:
        base_vcf, comp_vcf = params.base, params.comp
        regions["refined"] = original_stratify(base_vcf, comp_vcf, regions)
    else:
        base_vcf, comp_vcf = consolidate_bench_vcfs(args.benchdir)
        regions["refined"] = (regions["in_fn"] > 0) & (regions["in_fp"] > 0)
    logging.info("%d regions to be refined", regions["refined"].sum())

    reeval_bed = truvari.make_temp_filename(suffix=".bed")
    regions[regions["refined"]].to_csv(reeval_bed, sep='\t',
                                       header=False, index=False)

    # Send the vcfs to phab
    phab_vcf = os.path.join(args.benchdir, "phab.output.vcf.gz")
    to_eval_coords = (regions[regions["refined"]][["chrom", "start", "end"]]
                      .to_numpy()
                      .tolist())
    truvari.phab(to_eval_coords, base_vcf, args.reference, phab_vcf, buffer=PHAB_BUFFER,
                 mafft_params=args.mafft_params, comp_vcf=comp_vcf, prefix_comp=True,
                 threads=args.threads, method=args.align)

    # Now run bench on the phab harmonized variants
    logging.info("Running bench")
    matcher = truvari.Matcher(params)
    matcher.params.no_ref = 'a'
    outdir = os.path.join(args.benchdir, "phab_bench")
    m_bench = truvari.Bench(matcher, phab_vcf, phab_vcf, outdir, reeval_bed)
    m_bench.run()

    regions = refined_stratify(outdir, to_eval_coords, regions, args.threads)

    summary = (recount_variant_report(args.benchdir, outdir, regions)
               if args.recount else make_variant_report(regions))
    summary.clean_out()
    summary.write_json(os.path.join(
        args.benchdir, 'refine.variant_summary.json'))

    report = make_region_report(regions)
    regions.to_csv(os.path.join(
        args.benchdir, 'refine.regions.txt'), sep='\t', index=False)
    with open(os.path.join(args.benchdir, "refine.region_summary.json"), 'w') as fout:
        fout.write(json.dumps(report, indent=4))

    logging.info("Finished refine")
