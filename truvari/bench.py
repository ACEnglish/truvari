from __future__ import print_function
"""
Structural variant caller comparison tool
Given a benchmark and callset, calculate the recall/precision/f-measure
"""
import os
import re
import sys
import json
import bisect
import logging
import argparse

from truvari import *

from collections import defaultdict, namedtuple, OrderedDict

# External dependencies
import pysam
import pyfaidx

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
    thresg.add_argument("-P", "--pctsize", type=restricted_float, default=0.70,
                        help="Min pct allele size similarity (minvarsize/maxvarsize) (%(default)s)")
    thresg.add_argument("-O", "--pctovl", type=restricted_float, default=0.0,
                        help="Minimum pct reciprocal overlap (%(default)s)")
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

def annotate_tp(entry, score, pctsim, pctsize, pctovl, szdiff, stdist, endist, oentry):
    """
    Add the matching annotations to a vcf entry
    match_score, match_pctsim, match_pctsize, match_ovlpct, match_szdiff, \
                    match_stdist, match_endist, match_entry
    """
    entry.info["PctSeqSimilarity"] = pctsim
    entry.info["PctSizeSimilarity"] = pctsize
    entry.info["PctRecOverlap"] = pctovl
    entry.info["SizeDiff"] = szdiff
    entry.info["StartDistance"] = stdist
    entry.info["EndDistance"] = endist

def check_params(args):
    """
    Checks parameters as much as possible.
    All errors are written to stderr without logging since failures mean no output
    """
    check_fail = False
    if os.path.isdir(args.output):
        logging.error("Output directory '%s' already exists" % args.output)
        check_fail = True
    if not os.path.exists(args.comp):
        check_fail = True
        logging.error("File %s does not exist" % (args.comp))
    if not os.path.exists(args.base):
        check_fail = True
        logging.error("File %s does not exist" % (args.base))
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
    vcf_file = pysam.VariantFile(args.base)
    check_fail = False
    if sampleId is not None or sampleId not in vcf_file.header.samples:
        logging.error("Sample %s not found in vcf (%s)", sampleId, vcf_fn)
        check_fail = True
    if len(vcf_base.header.samples) == 0:
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
    
    check_fail = get_sample(args.base, args.bSample)
    check_fail = get_sample(args.comp, args.cSample)
    
    return check_fail

def setup_outputs(args):
    """
    Makes all of the output files
    return a ... to get to each of the
    """
    os.mkdir(args.output)
    setup_logging(args.debug, LogFileStderr(os.path.join(args.output, "log.txt")))
    logging.info("Params:\n%s", json.dumps(vars(args), indent=4))
    
    dat = namedtuple("outputs", ("vcf_base n_base_header sampleBase "
                                 "vcf_comp n_comp_header samplecomp "
                                 "tpb_out tpc_out fn_out fp_out "
                                 "b_filt c_filt stats_box"))

    vcf_base = pysam.VariantFile(args.base, 'r')
    n_base_header = edit_header(vcf_base)
    sampleBase = args.bSample if args.bSample else vcf_base.header.samples[0]
    
    vcf_comp = pysam.VariantFile(args.comp)
    n_comp_header = edit_header(vcf_comp)
    sampleComp = args.cSample if args.cSample else vcf_comp.header.samples[0]

    vcf_base = pysam.VariantFile(args.base)
    n_base_header = edit_header(vcf_base)
    # Setup outputs
    tpb_out = pysam.VariantFile(os.path.join(args.output, "tp-base.vcf"), 'w', header=n_base_header)
    tpc_out = pysam.VariantFile(os.path.join(args.output, "tp-call.vcf"), 'w', header=n_comp_header)

    fn_out = pysam.VariantFile(os.path.join(args.output, "fn.vcf"), 'w', header=n_base_header)
    fp_out = pysam.VariantFile(os.path.join(args.output, "fp.vcf"), 'w', header=n_comp_header)

    b_filt = pysam.VariantFile(os.path.join(args.output, "base-filter.vcf"), 'w', header=n_base_header)
    c_filt = pysam.VariantFile(os.path.join(args.output, "call-filter.vcf"), 'w', header=n_comp_header)

    stats_box = utils.StatsBox()
    
    return dat(vcf_base, None, sampleBase, vcf_comp, None, samplecomp, 
               tpb_out, tpc_out, fn_out, fp_out, b_filt, c_filt, stats_box)

def filter_call(entry, args, sizeA, outputs, base=True):
    """
    Given an entry, the parse_args arguments, and the entry's size
    check if it should be excluded from further analysis
    """
    prefix = "base " if base else "comp "
    if sizeA < args.sizemin or sizeA > args.sizemax:
        outputs.stats_box[prefix + "size filtered"] += 1
        if base:
            outputs.b_filt.write(entry)
        else:
            outputs.c_filt.write(entry)
        return True
    samp = outputs.sampleBase if base else outputs.sampleComp
    if args.no_ref in ["a", "b"] and not entry_is_variant(entry, samp):
        outputs.stats_box[prefix + "gt filtered"] += 1
        if base:
            outputs.b_filt.write(entry)
        else:
            outputs.c_filt.write(entry)
        return True
    if args.passonly and "PASS" not in entry.filter:
        if args.debug:
            logging.debug(which + "variant has no PASS FILTER and is being excluded from comparison - %s",
                              entry)
            return True
    return False

def write_fn(base_entry, outputs):
    """
    Write an entry to the false_negative output
    """
    # No overlaps, don't even bother checking
    n_base_entry = copy_entry(base_entry, outputs.n_base_header)
    n_base_entry.info["NumNeighbors"] = 0
    n_base_entry.info["NumThresholdNeighbors"] = 0
    outputs.stats_box["FN"] += 1
    outputs.fn_out.write(n_base_entry)

def match_calls(base_entry, comp_entry, astart, aend, sizeA, regions, args, outputs):
    """
    Compare the base and comp entries.
    We provied astart...sizeA because we've presumably calculated it before
    Note - This is the crucial component of matching.. so needs to be better 
    pulled apart for reusability and put into comparisons
    """
    sizeB = get_vcf_entry_size(comp_entry)
    if sizeB < args.sizefilt:
        return

    # Double ensure OVERLAP - there's a weird edge case where fetch with
    # the interval tree can return non-overlaps
    bstart, bend = get_vcf_boundaries(comp_entry)
    if not overlaps(astart - args.refdist, aend + args.refdist, bstart, bend):
        return
    
    if not regions.include(comp_entry):
        return
    
    # Someone in the Base call's neighborhood, we'll see if it passes comparisons
    num_neighbors += 1

    if args.no_ref in ["a", "c"] and not entry_is_variant(comp_entry, outputs.sampleBase):
        logging.debug("%s is uncalled", comp_entry)
        return

    if args.gtcomp and not gt_comp(base_entry, comp_entry, outputs.sampleBase, outputs.sampleComp):
        logging.debug("%s and %s are not the same genotype", str(base_entry), str(comp_entry))
        return

    if not args.typeignore and not same_variant_type(base_entry, comp_entry):
        logging.debug("%s and %s are not the same SVTYPE", str(base_entry), str(comp_entry))
        return

    size_similarity, size_diff = var_sizesim(sizeA, sizeB)
    if size_similarity < args.pctsize:
        logging.debug("%s and %s size similarity is too low (%f)", str(base_entry), 
                     str(comp_entry), size_similarity)
        return

    ovl_pct = get_rec_ovl(astart, aend, bstart, bend)
    if ovl_pct < args.pctovl:
        logging.debug("%s and %s overlap percent is too low (%f)", str(base_entry), str(comp_entry), ovl_pct)
        return

    if args.pctsim > 0:
        seq_similarity = var_pctsim_lev(base_entry, comp_entry, reference)
        if seq_similarity < args.pctsim:
            logging.debug("%s and %s sequence similarity is too low (%f)", str(
                base_entry), str(comp_entry), seq_similarity)
            return
    else:
        seq_similarity = 0

    start_distance = astart - bstart
    end_distance = aend - bend

    score = get_weighted_score(seq_similarity, size_similarity, ovl_pct)
    # If you put these numbers in an object, it'd be easier to pass round
    # You'd just need to make it sortable
    return (score, seq_similarity, size_similarity, ovl_pct, size_diff, start_distance, end_distance, comp_entry)

def output_base_match(base_entry, num_neighbors, thresh_neighbors, myid, matched_calls, outputs):
    """
    Writes a base call after it has gone through matching
    """
    base_entry = copy_entry(base_entry, outputs.n_base_header)
    base_entry.info["NumNeighbors"] = num_neighbors
    base_entry.info["NumThresholdNeighbors"] = len(thresh_neighbors)
    base_entry.info["MatchId"] = pbarcnt
    
    if len(thresh_neighbors) == 0:
        # False negative
        outputs.stats_box["FN"] += 1
        outputs.fn_out.write(base_entry)
        return

    logging.debug("Picking from candidate matches:\n%s", "\n".join([str(x) for x in thresh_neighbors]))
    truvari.match_sorter(thresh_neighbors)
    match_score, match_pctsim, match_pctsize, match_ovlpct, match_szdiff, \
        match_stdist, match_endist, match_entry = thresh_neighbors[0]
    logging.debug("Best match is %s", str(match_entry))

    base_entry.info["TruScore"] = match_score

    annotate_tp(base_entry, *thresh_neighbors[0])
    outputs.tpb_out.write(base_entry)

    # Don't double count calls found before
    b_key = vcf_to_key('b', base_entry)
    if not matched_calls[b_key]:
        # Interesting...
        outputs.stats_box["TP-base"] += 1
        if gt_comp(base_entry, match_entry, sampleBase, sampleComp):
            outputs.stats_box["TP-base_TP-gt"] += 1
        else:
            outputs.stats_box["TP-base_FP-gt"] += 1
    # Mark the call for multimatch checking
    matched_calls[b_key] = True
 
def report_best_match(base_entry, num_neighbors, thresh_neighbors, myid, matched_calls, outputs):
    """
    Pick and record the best base_entry
    """
    output_base_match(base_entry, num_neighbors, thresh_neighbors, myid, matched_calls, outputs)
    
    # Work through the comp calls
    for thresh_neighbor in thresh_neighbors:
        match_score, match_pctsim, match_pctsize, match_ovlpct, match_szdiff, \
            match_stdist, match_endist, match_entry = thresh_neighbor

        # Multimatch checking
        c_key = vcf_to_key('c', match_entry)
        if not matched_calls[c_key]: # unmatched
            outputs.stats_box["TP-call"] += 1
            if gt_comp(base_entry, match_entry, outputs.sampleBase, outputs.sampleComp):
                outputs.stats_box["TP-call_TP-gt"] += 1
            else:
                outputs.stats_box["TP-call_FP-gt"] += 1
        elif not args.multimatch:
            # Used this one and it can't multimatch
            continue
        
        logging.debug("Matching %s and %s", str(base_entry), str(match_entry))
        match_entry = copy_entry(match_entry, n_comp_header)
        match_entry.info["TruScore"] = match_score
        match_entry.info["NumNeighbors"] = num_neighbors
        match_entry.info["NumThresholdNeighbors"] = len(thresh_neighbors)
        match_entry.info["MatchId"] = pbarcnt
        annotate_tp(match_entry, *thresh_neighbor)
        tpc_out.write(match_entry)
        
        # Mark the call for multimatch checking
        matched_calls[c_key] = True
        if not args.multimatch:  # We're done here
            break

def parse_fps(matched_calls, tot_comp_entries, regions, args, outputs):
    """
    Report all the false-positives and comp filered calls
    """
    if args.prog:
        bar = setup_progressbar(tot_comp_entries)

    # Reset
    vcf_comp = pysam.VariantFile(args.comp)
    for cnt, entry in enumerate(regions.iterate(vcf_comp)):
        if args.passonly and 'PASS' not in entry.filter:
            continue
        if args.prog:
            bar.update(cnt + 1)
        if matched_calls[vcf_to_key('c', entry)]:
            continue
        size = get_vcf_entry_size(entry)
        if size < args.sizemin or size > args.sizemax:
            outputs.c_filt.write(copy_entry(entry, outputs.n_comp_header))
            outputs.stats_box["call size filtered"] += 1
        elif args.no_ref in ["a", "c"] and not entry_is_variant(entry, outputs.sampleB):
            outputs.stats_box["call gt filtered"] += 1
        elif regions.include(entry):
            outputs.fp_out.write(copy_entry(entry, outputs.n_comp_header))
            outputs.stats_box["FP"] += 1
    
    if args.prog:
        bar.finish()

def bench_main(cmdargs):
    """
    Main entry point for running Truvari Benchmarking
    """
    args = parse_args(cmdargs)
    
    if check_params(args) or check_iputs(args):
        sys.stderr("Couldn't run Truvari. Please fix parameters\n")
        exit(100)
    
    # We can now 'safely' perform everything
    outputs = setup_outputs()
    reference = pyfaidx.Fasta(args.reference) if args.reference else None
    
    logging.info("Creating call interval tree for overlap search")
    regions = GenomeTree(outputs.vcf_base, outputs.vcf_comp, args.includebed, args.sizemax)
    span_lookup, tot_comp_entries, cmp_entries = make_interval_tree(
        regions.iterate(outputs.vcf_comp), args.sizefilt, args.sizemax, args.passonly)
    logging.info("%d call variants in total", tot_comp_entries)
    logging.info("%d call variants within size range (%d, %d)", cmp_entries, args.sizefilt, args.sizemax)
    
    num_entries = 0
    bar = None
    if args.prog:
        for entry in regions.iterate(outputs.vcf_base):
            num_entries += 1
        logging.info("%s base variants", num_entries)
        bar = setup_progressbar(num_entries)
    
    # Reset
    outputs.vcf_base = pysam.VariantFile(args.base, 'r')
    outputs.n_base_header = edit_header(vcf_base)
    outputs.vcf_comp = pysam.VariantFile(args.comp) # not sure I have to do this one, though
    outputs.n_comp_header = edit_header(vcf_comp)

    # Calls that have been matched up
    matched_calls = defaultdict(bool)
    
    # for variant in base - do filtering on it and then try to match it to comp
    logging.info("Matching base to calls")

    for pbarcnt, base_entry in enumerate(regions.iterate(outputs.vcf_base)):
        if args.prog:
            bar.update(pbarcnt)
        
        sizeA = get_vcf_entry_size(base_entry)
        
        if filter_call(base_entry, args, sizeA):
            continue

        outputs.stats_box["base cnt"] += 1

        fetch_start, fetch_end = fetch_coords(span_lookup, base_entry, args.refdist)
        # No overlaps, don't even bother checking
        if fetch_start is None and fetch_end is None:
            write_fn(base_entry, outputs)
        
        # IntervalTree can give boundaries past REFDIST in the case of Inversions where start>end
        # We still need to fetch on the expanded boundaries so we can test them, but
        # we need to filter calls that otherwise shouldn't be considered
        # see the bstart/bend below
        astart, aend = get_vcf_boundaries(base_entry)

        # Search for comparison vcf entries as potential matches
        thresh_neighbors = []
        num_neighbors = 0
        
        # +- 1 just to be safe because why not
        for comp_entry in vcf_comp.fetch(base_entry.chrom, max(0, fetch_start - 1), fetch_end + 1):
            if filter_call(comp_entry, args, outputs, False):
                continue
                
            # There is a race condition here that could potentially mismatch things
            # If base1 passes matching call1 and then base2 passes matching call1
            # better, it can't use it and we mismatch 
            # UPDATE: by default we don't enforce one-match
            logging.debug("Comparing %s %s", str(base_entry), str(comp_entry))
            b_key = vcf_to_key('c', comp_entry)
            if not args.multimatch and matched_calls[b_key]:
                logging.debug("No match because comparison call already matched")
                continue
            
            dat = match_calls(base_entry, comp_entry, astart, aend, sizeA, regions, args, outputs)
            
            thresh_neighbors.append(dat)

        # Finished with this base entry
        report_best_match(base_entry, num_neighbors, thresh_neighbors, myid, matched_calls, outputs)
        
    if args.prog:
        bar.finish()
    
    # What is this... very broken is what
    outputs.stats_box.calc_performance(True)
    
    logging.info("Parsing FPs from calls")
    parse_fps(matched_calls, tot_comp_entries, regions, args, outputs)
    
    # call count is just of those used were used
    outputs.stats_box["call cnt"] = outputs.stats_box["TP-base"] + outputs.stats_box["FP"]

    # Close to flush vcfs
    outputs.tpb_out.close()
    outputs.b_filt.close()
    outputs.tpc_out.close()
    outputs.c_filt.close()
    outputs.fn_out.close()
    outputs.fp_out.close()
    
    # make stats
    with open(os.path.join(args.output, "summary.txt"), 'w') as fout:
        fout.write(json.dumps(outputs.stats_box, indent=4))
        logging.info("Stats: %s", json.dumps(outputs.stats_box, indent=4))

    if args.giabreport:
        make_giabreport(args, outputs.stats_box)

    logging.info("Finished")
