"""
Structural variant caller comparison tool
Given a benchmark and callset, calculate the recall/precision/f-measure
"""
import os
import sys
import copy
import json
import logging
import argparse
import itertools

from collections import defaultdict, OrderedDict, Counter

import pysam
import numpy as np

import truvari
from truvari.giab_report import make_giabreport

class StatsBox(OrderedDict):
    """
    Make a blank stats box for counting TP/FP/FN and calculating performance
    """

    def __init__(self):
        super().__init__()
        self["TP-base"] = 0
        self["TP-call"] = 0
        self["FP"] = 0
        self["FN"] = 0
        self["precision"] = 0
        self["recall"] = 0
        self["f1"] = 0
        self["base cnt"] = 0
        self["call cnt"] = 0
        self["TP-call_TP-gt"] = 0
        self["TP-call_FP-gt"] = 0
        self["TP-base_TP-gt"] = 0
        self["TP-base_FP-gt"] = 0
        self["gt_concordance"] = 0
        self["gt_matrix"] = defaultdict(Counter)

    def calc_performance(self):
        """
        Calculate the precision/recall
        """
        do_stats_math = True
        if self["TP-base"] == 0 and self["FN"] == 0:
            logging.warning("No TP or FN calls in base!")
            do_stats_math = False
        elif self["TP-call"] == 0 and self["FP"] == 0:
            logging.warning("No TP or FP calls in comp!")
            do_stats_math = False

        # Final calculations
        if do_stats_math:
            self["precision"] = float(self["TP-call"]) / \
                (self["TP-call"] + self["FP"])
            self["recall"] = float(self["TP-base"]) / \
                (self["TP-base"] + self["FN"])
            if self["TP-call_TP-gt"] + self["TP-call_FP-gt"] != 0:
                self["gt_concordance"] = float(self["TP-call_TP-gt"]) / (self["TP-call_TP-gt"] +
                                                                         self["TP-call_FP-gt"])

        # f-measure
        neum = self["recall"] * self["precision"]
        denom = self["recall"] + self["precision"]
        if denom != 0:
            self["f1"] = 2 * (neum / denom)

def compare_chunk(chunk):
    """
    Given a filtered chunk, (from chunker) compare all of the calls
    """
    matcher, chunk_dict, chunk_id = chunk
    logging.debug(f"Comparing chunk {chunk_dict}")

    # All FPs
    if len(chunk_dict['base']) == 0:
        fps = []
        for cid, c in enumerate(chunk_dict['comp']):
            ret = truvari.MatchResult()
            ret.comp = c
            ret.matid = f"{chunk_id}._.{cid}"
            fps.append(ret)
            logging.debug("All FP chunk -> {ret}")
        return fps

    # All FNs
    if len(chunk_dict['comp']) == 0:
        fns = []
        for bid, b in enumerate(chunk_dict['base']):
            ret = truvari.MatchResult()
            ret.base = b
            ret.matid = f"{chunk_id}.{bid}._"
            logging.debug(f"All FN chunk -> {ret}")
            fns.append(ret)
        return fns

    # Build all v all matrix
    # Should allow short-circuiting here.
    # If we find a perfect match, don't need to keep comparing.
    match_matrix = []
    for bid, b in enumerate(chunk_dict['base']):
        base_matches = []
        for cid, c in enumerate(chunk_dict['comp']):
            mat = matcher.build_match(b, c, f"{chunk_id}.{bid}.{cid}")
            logging.debug(f"Made mat -> {mat}")
            base_matches.append(mat)
        match_matrix.append(base_matches)

    # For building custom pickers
    match_matrix = np.array(match_matrix)
    if matcher.params.multimatch == 'KEEPALL':
        return match_matrix
    # Sort and annotate matches
    ret = []
    if matcher.params.multimatch:
        ret = pick_multi_matches(match_matrix)
    else:
        ret = pick_single_matches(match_matrix)
    return ret


def pick_multi_matches(match_matrix):
    """
    Given a numpy array of MatchResults
    Pick each base/comp call's best match
    """
    ret = []
    for b_max in match_matrix.max(axis=1):
        b_max = copy.copy(b_max)
        b_max.comp = None
        ret.append(b_max)

    for c_max in match_matrix.max(axis=0):
        c_max = copy.copy(c_max)
        c_max.base = None
        ret.append(c_max)
    return ret


def pick_single_matches(match_matrix):
    """
    Given a numpy array of MatchResults
    for each base, get its best match and record that the comp call can't be used anymore
    Once all best pairs are found, return the remaining
    """
    ret = []
    base_cnt, comp_cnt = match_matrix.shape
    match_matrix = np.ravel(match_matrix)
    match_matrix.sort()
    used_comp = set()
    used_base = set()
    for match in match_matrix[::-1]:
        # No more matches to find
        if base_cnt == 0 and comp_cnt == 0:
            break

        base_is_used = str(match.base) in used_base
        comp_is_used = str(match.comp) in used_comp
        # Only write the comp
        if base_cnt == 0 and not comp_is_used:
            to_process = copy.copy(match)
            to_process.base = None
            to_process.state = False
            to_process.multi = True
            comp_cnt -= 1
            used_comp.add(str(to_process.comp))
            ret.append(to_process)
        # Only write the base
        elif comp_cnt == 0 and not base_is_used:
            to_process = copy.copy(match)
            to_process.comp = None
            to_process.state = False
            to_process.multi = True
            base_cnt -= 1
            used_base.add(str(to_process.base))
            ret.append(to_process)
        # Write both
        elif not base_is_used and not comp_is_used:
            to_process = copy.copy(match)
            base_cnt -= 1
            used_base.add(str(to_process.base))
            comp_cnt -= 1
            used_comp.add(str(to_process.comp))
            ret.append(to_process)
    return ret


###############
# VCF editing #
###############
def edit_header(my_vcf):
    """
    Add INFO for new fields to vcf
    """
    header = my_vcf.header.copy()
    header.add_line(('##INFO=<ID=TruScore,Number=1,Type=Integer,'
                     'Description="Truvari score for similarity of match">'))
    header.add_line(('##INFO=<ID=PctSeqSimilarity,Number=1,Type=Float,'
                     'Description="Pct sequence similarity between this variant and its closest match">'))
    header.add_line(('##INFO=<ID=PctSizeSimilarity,Number=1,Type=Float,'
                     'Description="Pct size similarity between this variant and its closest match">'))
    header.add_line(('##INFO=<ID=PctRecOverlap,Number=1,Type=Float,'
                     'Description="Percent reciprocal overlap percent of the two calls\' coordinates">'))
    header.add_line(('##INFO=<ID=StartDistance,Number=1,Type=Integer,'
                     'Description="Distance of the base call\'s end from comparison call\'s start">'))
    header.add_line(('##INFO=<ID=EndDistance,Number=1,Type=Integer,'
                     'Description="Distance of the base call\'s end from comparison call\'s end">'))
    header.add_line(('##INFO=<ID=SizeDiff,Number=1,Type=Float,'
                     'Description="Difference in size of base and comp calls">'))
    header.add_line(('##INFO=<ID=GTMatch,Number=0,Type=Flag,'
                     'Description="Base/Comparison Genotypes match">'))
    header.add_line(('##INFO=<ID=MatchId,Number=1,Type=String,'
                     'Description="Id to help tie base/comp calls together {chunkid}.{baseid}.{compid}">'))
    header.add_line(('##INFO=<ID=Multi,Number=0,Type=Flag,'
                     'Description="Call is false due to non-multimatching">'))
    return header


def annotate_entry(entry, match, header):
    """
    Make a new entry with all the information
    """
    entry.translate(header)
    entry.info["PctSeqSimilarity"] = round(match.seqsim, 4) if match.seqsim is not None else None
    entry.info["PctSizeSimilarity"] = round(match.sizesim, 4) if match.sizesim is not None else None
    entry.info["PctRecOverlap"] = round(match.ovlpct, 4) if match.ovlpct is not None else None
    entry.info["SizeDiff"] = match.sizediff
    entry.info["StartDistance"] = match.st_dist
    entry.info["EndDistance"] = match.ed_dist
    entry.info["GTMatch"] = match.gt_match
    entry.info["TruScore"] = int(match.score) if match.score else None
    entry.info["MatchId"] = match.matid
    entry.info["Multi"] = match.multi


def output_writer(call, outs, sizemin):
    """
    Annotate a MatchResults' entries, write to the apppropriate file in outs
    and do the stats counting.
    Writer is responsible for handling FPs between sizefilt-sizemin
    """
    box = outs["stats_box"]
    if call.base:
        box["base cnt"] += 1
        annotate_entry(call.base, call, outs['n_base_header'])
        if call.state:
            gtBase = str(call.base_gt)
            gtComp = str(call.comp_gt)
            box["gt_matrix"][gtBase][gtComp] += 1

            box["TP-base"] += 1
            outs["tpb_out"].write(call.base)
            if call.gt_match:
                box["TP-base_TP-gt"] += 1
            else:
                box["TP-base_FP-gt"] += 1
        else:
            box["FN"] += 1
            outs["fn_out"].write(call.base)

    if call.comp:
        annotate_entry(call.comp, call, outs['n_comp_header'])
        if call.state:
            box["call cnt"] += 1
            box["TP-call"] += 1
            outs["tpc_out"].write(call.comp)
            if call.gt_match:
                box["TP-call_TP-gt"] += 1
            else:
                box["TP-call_FP-gt"] += 1
        elif truvari.entry_size(call.comp) >= sizemin:
            # The if is because we don't count FPs between sizefilt-sizemin
            box["call cnt"] += 1
            box["FP"] += 1
            outs["fp_out"].write(call.comp)


##########################
# Parameters and outputs #
##########################
def parse_args(args):
    """
    Pull the command line parameters
    """
    defaults = truvari.Matcher.make_match_params()
    parser = argparse.ArgumentParser(prog="bench", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--base", type=str, required=True,
                        help="Baseline truth-set calls")
    parser.add_argument("-c", "--comp", type=str, required=True,
                        help="Comparison set of calls")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output directory")
    parser.add_argument("-f", "--reference", type=str, default=None,
                        help="Indexed fasta used to call variants")
    parser.add_argument("--giabreport", action="store_true",
                        help="Parse output TPs/FNs for GIAB annotations and create a report")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Verbose logging")
    parser.add_argument("--prog", action="store_true",
                        help="Turn on progress monitoring")

    thresg = parser.add_argument_group("Comparison Threshold Arguments")
    thresg.add_argument("-r", "--refdist", type=truvari.restricted_int, default=defaults.refdist,
                        help="Max reference location distance (%(default)s)")
    thresg.add_argument("-u", "--unroll", action='store_true',
                        help="Use the unrolling procedure to perform sequence comparison")
    thresg.add_argument("-p", "--pctsim", type=truvari.restricted_float, default=defaults.pctsim,
                        help="Min percent allele sequence similarity. Set to 0 to ignore. (%(default)s)")
    thresg.add_argument("-B", "--minhaplen", type=truvari.restricted_int, default=defaults.minhaplen,
                        help="Minimum haplotype sequence length to create (%(default)s)")
    thresg.add_argument("-P", "--pctsize", type=truvari.restricted_float, default=defaults.pctsize,
                        help="Min pct allele size similarity (minvarsize/maxvarsize) (%(default)s)")
    thresg.add_argument("-O", "--pctovl", type=truvari.restricted_float, default=defaults.pctovl,
                        help="Min pct reciprocal overlap (%(default)s)")
    thresg.add_argument("-t", "--typeignore", action="store_true", default=defaults.typeignore,
                        help="Variant types don't need to match to compare (%(default)s)")
    thresg.add_argument("--dup-to-ins", action="store_true",
                        help="Assume DUP svtypes are INS (%(default)s)")
    thresg.add_argument("--use-lev", action="store_true",
                        help="Use the Levenshtein distance ratio instead of edlib editDistance ratio (%(default)s)")
    thresg.add_argument("-C", "--chunksize", type=truvari.restricted_int, default=defaults.chunksize,
                        help="Max reference distance to compare calls (%(default)s)")

    genoty = parser.add_argument_group("Genotype Comparison Arguments")
    genoty.add_argument("--gtcomp", action="store_true", default=defaults.gtcomp,
                        help="Compare GenoTypes of matching calls")
    genoty.add_argument("--bSample", type=str, default=None,
                        help="Baseline calls sample to use (first)")
    genoty.add_argument("--cSample", type=str, default=None,
                        help="Comparison calls sample to use (first)")

    filteg = parser.add_argument_group("Filtering Arguments")
    filteg.add_argument("-s", "--sizemin", type=truvari.restricted_int, default=defaults.sizemin,
                        help="Minimum variant size to consider from --comp (%(default)s)")
    filteg.add_argument("-S", "--sizefilt", type=truvari.restricted_int, default=None,
                        help="Minimum variant size to consider from --base (30)")
    filteg.add_argument("--sizemax", type=truvari.restricted_int, default=defaults.sizemax,
                        help="Maximum variant size to consider for comparison (%(default)s)")
    filteg.add_argument("--passonly", action="store_true", default=defaults.passonly,
                        help="Only consider calls with FILTER == PASS")
    filteg.add_argument("--no-ref", default=defaults.no_ref, choices=['a', 'b', 'c'],
                        help="Don't include 0/0 or ./. GT calls from all (a), base (b), or comp (c) vcfs (%(default)s)")
    filteg.add_argument("--includebed", type=str, default=None,
                        help="Bed file of regions in the genome to include only calls overlapping")
    filteg.add_argument("--extend", type=truvari.restricted_int, default=0,
                        help="Distance to allow comp entries outside of includebed regions (%(default)s)")
    filteg.add_argument("--multimatch", action="store_true", default=defaults.multimatch,
                        help=("Allow base calls to match multiple comparison calls, and vice versa. "
                              "Output vcfs will have redundant entries. (%(default)s)"))

    args = parser.parse_args(args)
    # When sizefilt is not provided and sizemin has been lowered below the default,
    # set sizefilt to sizemin. Otherwise, if sizefilt not provided, set to default
    # This just makes it easier to specify e.g. `-s 1` instead of `-s 1 -S 1`
    if args.sizefilt is None:
        if args.sizemin < defaults.sizefilt:
            args.sizefilt = args.sizemin
        else:
            args.sizefilt = defaults.sizefilt
    return args


def check_params(args):
    """
    Checks parameters as much as possible.
    All errors are written to stderr without logging since failures mean no output
    """
    check_fail = False
    if args.pctsim != 0 and not args.reference and not args.unroll:
        logging.error("--reference is required when --pctsim is set but not using --unroll")
        check_fail = True
    if args.chunksize < args.refdist:
        logging.error("--chunksize must be >= --refdist")
        check_fail = True
    if args.extend and args.includebed is None:
        logging.error("--extend can only be used when --includebed is set")
        check_fail = True
    if os.path.isdir(args.output):
        logging.error("Output directory '%s' already exists", args.output)
        check_fail = True
    if not os.path.exists(args.comp):
        logging.error("File %s does not exist", args.comp)
        check_fail = True
    if not os.path.exists(args.base):
        logging.error("File %s does not exist", args.base)
        check_fail = True
    if not args.comp.endswith(".gz"):
        logging.error(
            "Comparison vcf %s does not end with .gz. Must be bgzip'd", args.comp)
        check_fail = True
    if not os.path.exists(args.comp + '.tbi'):
        logging.error(
            "Comparison vcf index %s.tbi does not exist. Must be indexed", args.comp)
        check_fail = True
    if not args.base.endswith(".gz"):
        logging.error(
            "Base vcf %s does not end with .gz. Must be bgzip'd", args.base)
        check_fail = True
    if not os.path.exists(args.base + '.tbi'):
        logging.error(
            "Base vcf index %s.tbi does not exist. Must be indexed", args.base)
        check_fail = True
    if args.includebed and not os.path.exists(args.includebed):
        logging.error("Include bed %s does not exist", args.includebed)
        check_fail = True

    return check_fail


def check_sample(vcf_fn, sampleId=None):
    """
    Checks that a sample is inside a vcf
    Returns True if check failed
    """
    vcf_file = pysam.VariantFile(vcf_fn)
    check_fail = False
    if sampleId is not None and sampleId not in vcf_file.header.samples:
        logging.error("Sample %s not found in vcf (%s)", sampleId, vcf_fn)
        check_fail = True
    if len(vcf_file.header.samples) == 0:
        logging.error("No SAMPLE columns found in vcf (%s)", vcf_fn)
        check_fail = True
    return check_fail


def check_inputs(args):
    """
    Checks the inputs to ensure expected values are found inside of files
    Returns True if check failed
    """
    return check_sample(args.base, args.bSample) or check_sample(args.comp, args.cSample)


def setup_outputs(args):
    """
    Makes all of the output files
    Places the data into the shared space
    """
    os.mkdir(args.output)
    truvari.setup_logging(args.debug, truvari.LogFileStderr(
        os.path.join(args.output, "log.txt")))
    logging.info("Params:\n%s", json.dumps(vars(args), indent=4))
    logging.info(f"Truvari version: {truvari.__version__}")

    outputs = {}
    outputs["vcf_base"] = pysam.VariantFile(args.base)
    outputs["n_base_header"] = edit_header(outputs["vcf_base"])

    outputs["vcf_comp"] = pysam.VariantFile(args.comp)
    outputs["n_comp_header"] = edit_header(outputs["vcf_comp"])

    outputs["tpb_out"] = pysam.VariantFile(os.path.join(
        args.output, "tp-base.vcf"), 'w', header=outputs["n_base_header"])
    outputs["tpc_out"] = pysam.VariantFile(os.path.join(
        args.output, "tp-call.vcf"), 'w', header=outputs["n_comp_header"])

    outputs["fn_out"] = pysam.VariantFile(os.path.join(
        args.output, "fn.vcf"), 'w', header=outputs["n_base_header"])
    outputs["fp_out"] = pysam.VariantFile(os.path.join(
        args.output, "fp.vcf"), 'w', header=outputs["n_comp_header"])

    outputs["stats_box"] = StatsBox()
    return outputs


def close_outputs(outputs):
    """
    Close all the files
    """
    outputs["tpb_out"].close()
    outputs["tpc_out"].close()
    outputs["fn_out"].close()
    outputs["fp_out"].close()


def bench_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)

    if check_params(args) or check_inputs(args):
        sys.stderr.write("Couldn't run Truvari. Please fix parameters\n")
        sys.exit(100)

    matcher = truvari.Matcher(args=args)
    outputs = setup_outputs(args)

    base = pysam.VariantFile(args.base)
    comp = pysam.VariantFile(args.comp)

    regions = truvari.RegionVCFIterator(base, comp,
                                        args.includebed,
                                        args.sizemax)

    regions.merge_overlaps()
    regions_extended = regions.extend(args.extend) if args.extend else regions

    base_i = regions.iterate(base)
    comp_i = regions_extended.iterate(comp)

    chunks = truvari.chunker(matcher, ('base', base_i), ('comp', comp_i))
    for call in itertools.chain.from_iterable(map(compare_chunk, chunks)):
        # setting non-matched call variants that are not fully contained in the original regions to None
        # These don't count as FP or TP and don't appear in the output vcf files
        if args.extend and call.comp is not None and not call.state and not regions.include(call.comp):
            call.comp = None
        output_writer(call, outputs, args.sizemin)

    #with concurrent.futures.ThreadPoolExecutor(max_workers = args.threads) as executor:
    #    future_chunks = {executor.submit(compare_chunk, c): c for c in chunks}
    #    for future in concurrent.futures.as_completed(future_chunks):
    #        result = None
    #        try:
    #            result = future.result()
    #        except Exception as exc:
    #            logging.error('%r generated an exception: %s', result, exc)
    #        else:
    #            for call in result: #future_chunks[future]:

    with open(os.path.join(args.output, "summary.txt"), 'w') as fout:
        box = outputs["stats_box"]
        box.calc_performance()
        fout.write(json.dumps(box, indent=4))
        logging.info("Stats: %s", json.dumps(box, indent=4))

    close_outputs(outputs)

    if args.giabreport:
        make_giabreport(args, outputs["stats_box"])

    logging.info("Finished bench")
