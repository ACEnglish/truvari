"""
Structural variant caller comparison tool
Given a benchmark and callset, calculate the recall/precision/f-measure
"""
# pylint: disable=too-many-statements, no-member
import os
import sys
import json
import logging
import argparse
from collections import defaultdict, namedtuple

import pysam
import pyfaidx

import truvari

MATCHRESULT = namedtuple("matchresult", ("score seq_similarity size_similarity "
                                         "ovl_pct size_diff start_distance "
                                         "end_distance match_entry"))

def parse_args(args):
    """
    Pull the command line parameters
    """
    def restricted_float(x):
        x = float(x)
        if x < 0.0 or x > 1.0:
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
        return x

    parser = argparse.ArgumentParser(prog="truvari", description=__doc__,
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
    thresg.add_argument("-r", "--refdist", type=int, default=500,
                        help="Max reference location distance (%(default)s)")
    thresg.add_argument("-p", "--pctsim", type=restricted_float, default=0.70,
                        help="Min percent allele sequence similarity. Set to 0 to ignore. (%(default)s)")
    thresg.add_argument("-B", "--buffer", type=restricted_float, default=0.10,
                        help="Percent of the reference span to buffer the haplotype sequence created")
    thresg.add_argument("-P", "--pctsize", type=restricted_float, default=0.70,
                        help="Min pct allele size similarity (minvarsize/maxvarsize) (%(default)s)")
    thresg.add_argument("-O", "--pctovl", type=restricted_float, default=0.0,
                        help="Minimum pct reciprocal overlap (%(default)s) for DEL events")
    thresg.add_argument("-t", "--typeignore", action="store_true", default=False,
                        help="Variant types don't need to match to compare (%(default)s)")

    genoty = parser.add_argument_group("Genotype Comparison Arguments")
    genoty.add_argument("--gtcomp", action="store_true", default=False,
                        help="Compare GenoTypes of matching calls")
    genoty.add_argument("--bSample", type=str, default=None,
                        help="Baseline calls sample to use (first)")
    genoty.add_argument("--cSample", type=str, default=None,
                        help="Comparison calls sample to use (first)")

    filteg = parser.add_argument_group("Filtering Arguments")
    filteg.add_argument("-s", "--sizemin", type=int, default=50,
                        help="Minimum variant size to consider for comparison (%(default)s)")
    filteg.add_argument("-S", "--sizefilt", type=int, default=30,
                        help="Minimum variant size to load into IntervalTree (%(default)s)")
    filteg.add_argument("--sizemax", type=int, default=50000,
                        help="Maximum variant size to consider for comparison (%(default)s)")
    filteg.add_argument("--passonly", action="store_true", default=False,
                        help="Only consider calls with FILTER == PASS")
    filteg.add_argument("--no-ref", default=False, choices=['a', 'b', 'c'],
                        help="Don't include 0/0 or ./. GT calls from all (a), base (b), or comp (c) vcfs (%(default)s)")
    filteg.add_argument("--includebed", type=str, default=None,
                        help="Bed file of regions in the genome to include only calls overlapping")
    filteg.add_argument("--multimatch", action="store_true", default=False,
                        help=("Allow base calls to match multiple comparison calls, and vice versa. "
                              "Output vcfs will have redundant entries. (%(default)s)"))

    args = parser.parse_args(args)
    if args.pctsim != 0 and not args.reference:
        parser.error("--reference is required when --pctsim is set")

    return args

def edit_header(my_vcf):
    """
    Add INFO for new fields to vcf
    #Probably want to put in the PG whatever, too
    """
    # Update header
    # Edit Header
    header = my_vcf.header.copy()
    header.add_line(('##INFO=<ID=TruScore,Number=1,Type=Float,'
                     'Description="Truvari score for similarity of match">'))
    header.add_line(('##INFO=<ID=PctSeqSimilarity,Number=1,Type=Float,'
                     'Description="Pct sequence similarity between this variant and its closest match">'))
    header.add_line(('##INFO=<ID=PctSizeSimilarity,Number=1,Type=Float,'
                     'Description="Pct size similarity between this variant and its closest match">'))
    header.add_line(('##INFO=<ID=PctRecOverlap,Number=1,Type=Float,'
                     'Description="Percent reciprocal overlap percent of the two calls\' coordinates">'))
    header.add_line(('##INFO=<ID=StartDistance,Number=1,Type=Integer,'
                     'Description="Distance of this call\'s start from comparison call\'s start">'))
    header.add_line(('##INFO=<ID=EndDistance,Number=1,Type=Integer,'
                     'Description="Distance of this call\'s start from comparison call\'s start">'))
    header.add_line(('##INFO=<ID=SizeDiff,Number=1,Type=Float,'
                     'Description="Difference in size(basecall) and size(evalcall)">'))
    header.add_line(('##INFO=<ID=NumNeighbors,Number=1,Type=Integer,'
                     'Description="Number of calls in B that were in the neighborhood (REFDIST) of this call">'))
    header.add_line(('##INFO=<ID=NumThresholdNeighbors,Number=1,Type=Integer,'
                     'Description="Number of calls in B that are within threshold distances of this call">'))
    header.add_line(('##INFO=<ID=MatchId,Number=1,Type=Integer,'
                     'Description="Truvari uid to help tie tp-base.vcf and tp-call.vcf entries together">'))
    return header

def annotate_tp(entry, match_result):
    """
    Add the matching annotations to a vcf entry
    match_score, match_pctsim, match_pctsize, match_ovlpct, match_szdiff, \
                    match_stdist, match_endist, match_entry
    """
    entry.info["PctSeqSimilarity"] = match_result.seq_similarity
    entry.info["PctSizeSimilarity"] = match_result.size_similarity
    entry.info["PctRecOverlap"] = match_result.ovl_pct
    entry.info["SizeDiff"] = match_result.size_diff
    entry.info["StartDistance"] = match_result.start_distance
    entry.info["EndDistance"] = match_result.end_distance

def check_params(args):
    """
    Checks parameters as much as possible.
    All errors are written to stderr without logging since failures mean no output
    """
    check_fail = False
    if os.path.isdir(args.output):
        logging.error("Output directory '%s' already exists", args.output)
        check_fail = True
    if not os.path.exists(args.comp):
        check_fail = True
        logging.error("File %s does not exist", args.comp)
    if not os.path.exists(args.base):
        check_fail = True
        logging.error("File %s does not exist", args.base)
    if not args.comp.endswith(".gz"):
        check_fail = True
        logging.error("Comparison vcf %s does not end with .gz. Must be bgzip'd", args.comp)
    if not os.path.exists(args.comp + '.tbi'):
        check_fail = True
        logging.error("Comparison vcf index %s.tbi does not exist. Must be indexed", args.comp)

    return check_fail

def check_sample(vcf_fn, sampleId=None):
    """
    Return the given sampleId from the var_file
    if sampleId is None, return the first sample
    if there is no first sample in the var_file, raise an error
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
    Checks the inputs against the arguments as much as possible before creating any output
    Returns:
        vcf_bse
    """
    check_fail = False
    check_fail = check_sample(args.base, args.bSample)
    check_fail = check_sample(args.comp, args.cSample)
    return check_fail

def setup_outputs(args):
    """
    Makes all of the output files
    return a ... to get to each of the
    """
    os.mkdir(args.output)
    truvari.setup_logging(args.debug, truvari.LogFileStderr(os.path.join(args.output, "log.txt")))
    logging.info("Params:\n%s", json.dumps(vars(args), indent=4))

    outputs = {}

    outputs["vcf_base"] = pysam.VariantFile(args.base, 'r')
    outputs["n_base_header"] = edit_header(outputs["vcf_base"])
    outputs["sampleBase"] = args.bSample if args.bSample else outputs["vcf_base"].header.samples[0]

    outputs["vcf_comp"] = pysam.VariantFile(args.comp)
    outputs["n_comp_header"] = edit_header(outputs["vcf_comp"])
    outputs["sampleComp"] = args.cSample if args.cSample else outputs["vcf_comp"].header.samples[0]

    outputs["vcf_base"] = pysam.VariantFile(args.base)
    outputs["n_base_header"] = edit_header(outputs["vcf_base"])
    # Setup outputs
    outputs["tpb_out"] = pysam.VariantFile(os.path.join(args.output, "tp-base.vcf"), 'w', header=outputs["n_base_header"])
    outputs["tpc_out"] = pysam.VariantFile(os.path.join(args.output, "tp-call.vcf"), 'w', header=outputs["n_comp_header"])

    outputs["fn_out"] = pysam.VariantFile(os.path.join(args.output, "fn.vcf"), 'w', header=outputs["n_base_header"])
    outputs["fp_out"] = pysam.VariantFile(os.path.join(args.output, "fp.vcf"), 'w', header=outputs["n_comp_header"])

    outputs["stats_box"] = truvari.StatsBox()

    return outputs

def filter_call(entry, sizeA, sizemin, sizemax, no_ref, passonly, outputs, base=True):
    """
    Given an entry, the parse_args arguments, and the entry's size
    check if it should be excluded from further analysis
    """
    prefix = "base " if base else "call "
    if sizeA < sizemin or sizeA > sizemax:
        return True
    
    samp = outputs["sampleBase"] if base else outputs["sampleComp"]
    if no_ref in ["a", "b"] and not truvari.entry_is_variant(entry, samp):
        return True

    if passonly and truvari.filter_value(entry):
        logging.debug("%s variant has no PASS FILTER and is being excluded from comparison - %s",
                      prefix, entry)
        return True
    return False

def write_fn(base_entry, outputs):
    """
    Write an entry to the false_negative output
    """
    # No overlaps, don't even bother checking
    n_base_entry = truvari.copy_entry(base_entry, outputs["n_base_header"])
    n_base_entry.info["NumNeighbors"] = 0
    n_base_entry.info["NumThresholdNeighbors"] = 0
    outputs["stats_box"]["FN"] += 1
    outputs["fn_out"].write(n_base_entry)

def match_calls(base_entry, comp_entry, astart, aend, sizeA, sizeB, regions, reference, args, outputs):
    """
    Compare the base and comp entries.
    We provied astart...sizeA because we've presumably calculated it before
    Note - This is the crucial component of matching.. so needs to be better
    pulled apart for reusability and put into comparisons
    """
    if sizeB < args.sizefilt:
        return False

    # Double ensure OVERLAP - there's a weird edge case where fetch with
    # the interval tree can return non-overlaps
    bstart, bend = truvari.entry_boundaries(comp_entry)
    if not truvari.overlaps(astart - args.refdist, aend + args.refdist, bstart, bend):
        return False

    if not regions.include(comp_entry):
        return False

    # Someone in the Base call's neighborhood, we'll see if it passes comparisons

    if args.no_ref in ["a", "c"] and not truvari.entry_is_variant(comp_entry, outputs["sampleBase"]):
        logging.debug("%s is uncalled", comp_entry)
        return True

    if args.gtcomp and not truvari.entry_gt_comp(base_entry, comp_entry, outputs["sampleBase"], outputs["sampleComp"]):
        logging.debug("%s and %s are not the same genotype", str(base_entry), str(comp_entry))
        return True

    if not args.typeignore and not truvari.same_variant_type(base_entry, comp_entry):
        logging.debug("%s and %s are not the same SVTYPE", str(base_entry), str(comp_entry))
        return True

    size_similarity, size_diff = truvari.sizesim(sizeA, sizeB)
    if size_similarity < args.pctsize:
        logging.debug("%s and %s size similarity is too low (%f)", str(base_entry),
                      str(comp_entry), size_similarity)
        return True

    ovl_pct = truvari.reciprocal_overlap(astart, aend, bstart, bend)
    if truvari.entry_variant_type(base_entry) == "DEL" and ovl_pct < args.pctovl:
        logging.debug("%s and %s overlap percent is too low (%f)", str(base_entry), str(comp_entry), ovl_pct)
        return True

    if args.pctsim > 0:
        seq_similarity = truvari.entry_pctsim_lev(base_entry, comp_entry, reference, buf_len=args.buffer)
        if seq_similarity < args.pctsim:
            logging.debug("%s and %s sequence similarity is too low (%f)", str(
                base_entry), str(comp_entry), seq_similarity)
            return True
    else:
        seq_similarity = 0

    start_distance = astart - bstart
    end_distance = aend - bend

    score = truvari.weighted_score(seq_similarity, size_similarity, ovl_pct)

    return MATCHRESULT(score, seq_similarity, size_similarity, ovl_pct, size_diff,
                       start_distance, end_distance, comp_entry)

def output_base_match(base_entry, num_neighbors, thresh_neighbors, myid, matched_calls, outputs):
    """
    Writes a base call after it has gone through matching
    """
    base_entry = truvari.copy_entry(base_entry, outputs["n_base_header"])
    base_entry.info["NumNeighbors"] = num_neighbors
    base_entry.info["NumThresholdNeighbors"] = len(thresh_neighbors)
    base_entry.info["MatchId"] = myid

    if len(thresh_neighbors) == 0:
        # False negative
        outputs["stats_box"]["FN"] += 1
        outputs["fn_out"].write(base_entry)
        return

    logging.debug("Picking from candidate matches:\n%s", "\n".join([str(x) for x in thresh_neighbors]))
    truvari.match_sorter(thresh_neighbors)
    logging.debug("Best match is %s", str(thresh_neighbors[0].score))
    base_entry.info["TruScore"] = thresh_neighbors[0].score

    annotate_tp(base_entry, thresh_neighbors[0])
    outputs["tpb_out"].write(base_entry)

    # Don't double count calls found before
    b_key = truvari.entry_to_key('b', base_entry)
    if not matched_calls[b_key]:
        # Interesting...
        outputs["stats_box"]["TP-base"] += 1
        if truvari.entry_gt_comp(base_entry, thresh_neighbors[0].match_entry, outputs["sampleBase"], outputs["sampleComp"]):
            outputs["stats_box"]["TP-base_TP-gt"] += 1
        else:
            outputs["stats_box"]["TP-base_FP-gt"] += 1
    # Mark the call for multimatch checking
    matched_calls[b_key] = True

def report_best_match(base_entry, num_neighbors, thresh_neighbors, myid, matched_calls, outputs, args):
    """
    Pick and record the best base_entry
    """
    output_base_match(base_entry, num_neighbors, thresh_neighbors, myid, matched_calls, outputs)

    # Work through the comp calls
    for neigh in thresh_neighbors:
        # Multimatch checking
        c_key = truvari.entry_to_key('c', neigh.match_entry)
        if not matched_calls[c_key]: # unmatched
            outputs["stats_box"]["TP-call"] += 1
            if truvari.entry_gt_comp(base_entry, neigh.match_entry, outputs["sampleBase"], outputs["sampleComp"]):
                outputs["stats_box"]["TP-call_TP-gt"] += 1
            else:
                outputs["stats_box"]["TP-call_FP-gt"] += 1
        elif not args.multimatch:
            # Used this one and it can't multimatch
            continue

        logging.debug("Matching %s and %s", str(base_entry), str(neigh.match_entry))
        match_entry = truvari.copy_entry(neigh.match_entry, outputs["n_comp_header"])
        match_entry.info["TruScore"] = neigh.score
        match_entry.info["NumNeighbors"] = num_neighbors
        match_entry.info["NumThresholdNeighbors"] = len(thresh_neighbors)
        match_entry.info["MatchId"] = myid
        annotate_tp(match_entry, neigh)
        outputs["tpc_out"].write(match_entry)

        # Mark the call for multimatch checking
        matched_calls[c_key] = True
        if not args.multimatch:  # We're done here
            break

def parse_fps(matched_calls, tot_comp_entries, regions, args, outputs):
    """
    Report all the false-positives and comp filered calls
    """
    if args.prog:
        pbar = truvari.setup_progressbar(tot_comp_entries)

    # Reset
    vcf_comp = pysam.VariantFile(args.comp)
    for cnt, entry in enumerate(regions.iterate(vcf_comp)):
        # Here
        if args.prog:
            pbar.update(cnt + 1)
        
        size = truvari.entry_size(entry)
        if filter_call(entry, size, args.sizemin, args.sizemax, args.no_ref, args.passonly, outputs, False):
            continue
        
        if matched_calls[truvari.entry_to_key('c', entry)]:
            continue
            
        if regions.include(entry):
            outputs["fp_out"].write(truvari.copy_entry(entry, outputs["n_comp_header"]))
            outputs["stats_box"]["FP"] += 1

    if args.prog:
        pbar.finish()

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
    Main entry point for running Truvari Benchmarking
    """
    args = parse_args(cmdargs)

    if check_params(args) or check_inputs(args):
        sys.stderr.write("Couldn't run Truvari. Please fix parameters\n")
        sys.exit(100)

    # We can now 'safely' perform everything
    outputs = setup_outputs(args)
    reference = pysam.FastaFile(args.reference) if args.reference else None

    logging.info("Creating call interval tree for overlap search")
    regions = truvari.GenomeTree(outputs["vcf_base"], outputs["vcf_comp"], args.includebed, args.sizemax)
    span_lookup, tot_comp_entries, cmp_entries = truvari.make_interval_tree(
                regions.iterate(outputs["vcf_comp"]), args.sizefilt, args.sizemax, args.passonly)
    logging.info("%d call variants in total", tot_comp_entries)
    logging.info("%d call variants within size range (%d, %d)", cmp_entries, args.sizefilt, args.sizemax)

    num_entries = 0
    pbar = None
    if args.prog:
        for _ in regions.iterate(outputs["vcf_base"]):
            num_entries += 1
        logging.info("%s base variants", num_entries)
        pbar = truvari.setup_progressbar(num_entries)

    # Reset
    outputs["vcf_base"] = pysam.VariantFile(args.base, 'r')
    outputs["n_base_header"] = edit_header(outputs["vcf_base"])
    outputs["vcf_comp"] = pysam.VariantFile(args.comp)
    outputs["n_comp_header"] = edit_header(outputs["vcf_comp"])

    # Calls that have been matched up
    matched_calls = defaultdict(bool)

    # for variant in base - do filtering on it and then try to match it to comp
    logging.info("Matching base to calls")
    for pbarcnt, base_entry in enumerate(regions.iterate(outputs["vcf_base"])):
        if args.prog:
            pbar.update(pbarcnt)

        sizeA = truvari.entry_size(base_entry)

        if filter_call(base_entry, sizeA, args.sizemin, args.sizemax, args.no_ref, args.passonly, outputs, True):
            continue

        outputs["stats_box"]["base cnt"] += 1

        fetch_start, fetch_end = truvari.fetch_coords(span_lookup, base_entry, args.refdist)
        # No overlaps, don't even bother checking
        if fetch_start is None and fetch_end is None:
            write_fn(base_entry, outputs)
            continue

        # IntervalTree can give boundaries past REFDIST in the case of Inversions where start>end
        # We still need to fetch on the expanded boundaries so we can test them, but
        # we need to filter calls that otherwise shouldn't be considered
        # see the bstart/bend below
        astart, aend = truvari.entry_boundaries(base_entry)

        # Search for comparison vcf entries as potential matches
        thresh_neighbors = []
        num_neighbors = 0

        # +- 1 just to be safe because why not
        for comp_entry in outputs["vcf_comp"].fetch(base_entry.chrom, max(0, fetch_start - 1), fetch_end + 1):
            sizeB = truvari.entry_size(comp_entry)
            if filter_call(comp_entry, sizeB, args.sizefilt, args.sizemax, args.no_ref, args.passonly, outputs, False):
                continue

            # There is a race condition here that could potentially mismatch things
            # If base1 passes matching call1 and then base2 passes matching call1
            # better, it can't use it and we mismatch
            # UPDATE: by default we don't enforce one-match
            logging.debug("Comparing %s %s", str(base_entry), str(comp_entry))
            if not args.multimatch and matched_calls[truvari.entry_to_key('c', comp_entry)]:
                logging.debug("No match because comparison call already matched")
                continue
            mat = match_calls(base_entry, comp_entry, astart, aend, sizeA, sizeB, regions,
                              reference, args, outputs)
            if mat:
                num_neighbors += 1
            else:
                continue
            if isinstance(mat, bool):
                continue
            thresh_neighbors.append(mat)

        # Finished with this base entry
        report_best_match(base_entry, num_neighbors, thresh_neighbors, pbarcnt, matched_calls, outputs, args)

    if args.prog:
        pbar.finish()

    outputs["stats_box"].calc_performance(True)

    logging.info("Parsing FPs from calls")
    parse_fps(matched_calls, tot_comp_entries, regions, args, outputs)

    # call count is just of those used were used
    outputs["stats_box"]["call cnt"] = outputs["stats_box"]["TP-base"] + outputs["stats_box"]["FP"]

    # Close to flush vcfs
    close_outputs(outputs)

    # make stats
    outputs["stats_box"].calc_performance(False)
    with open(os.path.join(args.output, "summary.txt"), 'w') as fout:
        fout.write(json.dumps(outputs["stats_box"], indent=4))
        logging.info("Stats: %s", json.dumps(outputs["stats_box"], indent=4))

    if args.giabreport:
        truvari.make_giabreport(args, outputs["stats_box"])

    logging.info("Finished")
