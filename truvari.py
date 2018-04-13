#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import json
import bisect
import logging
import argparse
import warnings

from collections import defaultdict
# External dependencies

import vcf
import swalign
import Levenshtein
import progressbar
from intervaltree import IntervalTree

USAGE = """Structural variant caller comparison tool
Given a benchmark and callset, calculate the recall/precision/f-measure"""


# setup logging - need to write to file AND stderr
class LogFileStderr(object):

    """ Write to stderr and a file"""

    def __init__(self, fn):
        """ keep these props """
        self.name = fn
        self.file_handler = open(fn, 'w')

    def write(self, *args):
        """ Write to both """
        sys.stderr.write(*args)
        self.file_handler.write(*args)

    def flush(self):
        """ Flush both """
        sys.stderr.flush()
        self.file_handler.flush()


def setup_logging(debug=False, stream=sys.stderr, log_format=None):
    """
    Default logger
    """
    logLevel = logging.DEBUG if debug else logging.INFO
    if log_format is None:
        log_format = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(stream=stream, level=logLevel, format=log_format)
    logging.info("Running %s", " ".join(sys.argv))

    def sendWarningsToLog(message, category, filename, lineno):
        """
        Put warnings into logger
        """
        logging.warning('%s:%s: %s:%s', filename, lineno, category.__name__, message)
        return
    # pylint: disable=unused-variable
    old_showwarning = warnings.showwarning
    warnings.showwarning = sendWarningsToLog


def make_interval_tree(vcf_file, sizemin=10, passonly=False):
    """
    Return a dictonary of {chr:[start, end,..], ...}
    that can be queried with bisect to return a span to fetch variants against
    could stand to make a progress bar since this now takes so long
    """
    n_entries = 0
    cmp_entries = 0
    lookup = defaultdict(IntervalTree)
    for entry in vcf_file:
        if passonly and len(entry.FILTER):
            continue
        n_entries += 1
        start, end = get_vcf_boundaries(entry)
        if get_vcf_entry_size(entry) < sizemin:
            continue
        cmp_entries += 1
        lookup[entry.CHROM].addi(start, end, entry.start)
    return lookup, n_entries, cmp_entries


def vcf_to_key(entry):
    """
    Turn a vcf entry into a hashable key
    chr:pos:ref:alt
    helpful for not re-using variants
    BUG: if a caller redundantly calls a variant. It will be collapsed
    """
    start, end = get_vcf_boundaries(entry)
    return "%s:%d-%d(%s|%s)" % (entry.CHROM, start, end, entry.REF, str(entry.ALT[0]))


def var_sizesim(entryA, entryB, type_match=True):
    """
    Calculate the size similarity pct for the two entries
    compares the longer of entryA's two alleles (REF or ALT)
    """
    #Again, needs work
    if type_match and (len(entryA.REF) < len(entryA.ALT[0])) != (len(entryB.REF) < len(entryB.ALT[0])):
        return 0.0

    sizeA = get_vcf_entry_size(entryA)
    sizeB = get_vcf_entry_size(entryB)
    return min(sizeA, sizeB)/float(max(sizeA, sizeB))

def do_swalign(seq1, seq2, match=2, mismatch=-1, gap_penalty=-2, gap_extension_decay=0.5):
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring, gap_penalty=gap_penalty, gap_extension_decay=gap_extension_decay)
    aln = sw.align(seq1, seq2)
    return aln

def var_pctsim_lev(entryA, entryB, type_match=True):
    """
    Use Levenshtein distance ratio as the pct similarity
    """
    #This is broken... 
    if type_match and (len(entryA.REF) < len(entryA.ALT)) != (len(entryB.REF) < len(entryB.ALT)):
        return 0.0
    # Shortcut to save compute - probably unneeded
    if entryA.REF == entryB.REF and entryA.ALT[0] == entryB.ALT[0]:
        return 1.0
    
    if len(entryA.REF) < len(entryA.ALT[0]):
        return Levenshtein.ratio(str(entryA.ALT[0]), str(entryB.ALT[0]))
    return Levenshtein.ratio(str(entryA.REF), str(entryB.REF))
    

def var_pctsim_sw(entryA, entryB, type_match=True):
    """
    Find all the entries in vcfB that are within max_dist of entryA
    if type_match, only consider variants of the same type
    Works using swalign instead of edlib
    Returns the percent similarity
    """
    #This is broken... 
    if type_match and (len(entryA.REF) < len(entryA.ALT)) != (len(entryB.REF) < len(entryB.ALT)):
        return 0.0
    # Shortcut to save compute
    if entryA.REF == entryB.REF and entryA.ALT[0] == entryB.ALT[0]:
        return 1.0
    ref_aln = do_swalign(str(entryA.REF), str(entryB.REF))
    alt_aln = do_swalign(str(entryA.ALT[0]), str(entryB.ALT[0]))
    mat_tot = ref_aln.matches + alt_aln.matches
    mis_tot = ref_aln.mismatches + ref_aln.mismatches
    denom = float(mis_tot + mat_tot)
    if denom == 0:
        return 0
    ident = mat_tot / denom
    return ident

def __defunct_var_pctsim_ed(entryA, entryB, type_match=True):
    """
    Find all the entries in vcfB that are within max_dist of entryA
    if type_match, only consider variants of the same type

    Returns the percent similarity
    """
    #This is broken... 
    if type_match and (len(entryA.REF) < len(entryA.ALT)) != (len(entryB.REF) < len(entryB.ALT)):
        return 0.0
    # Shortcut to save compute
    if entryA.REF == entryB.REF and entryA.ALT[0] == entryB.ALT[0]:
        return 1.0
    ref_dist = edlib.align(str(entryA.REF), str(entryB.REF))["editDistance"]
    max_ref = max(len(entryA.REF), len(entryB.REF))
    # Only consider the first allele. This is a problem
    alt_dist = edlib.align(str(entryA.ALT[0]), str(entryB.ALT[0]))["editDistance"]
    max_alt = max(len(str(entryA.ALT[0])), len(str(entryB.ALT[0])))

    # dumbly, we'll just calculate some length and similarity percent
    ed = 1 - ((ref_dist + alt_dist) / float(max_ref + max_alt))
    return ed

def overlaps(s1, e1, s2, e2):
    """
    Check if two ranges have overlap
    """
    s_cand = max(s1, s2)
    e_cand = min(e1, e2)
    return s_cand < e_cand

def fetch_coords(lookup, entry, dist=0):
    """
    Get the minimum/maximum fetch coordinates to find all variants within dist of variant
    """

    start, end = get_vcf_boundaries(entry)
    start -= dist
    end += dist
    # Membership queries are fastest O(1)
    if not lookup[entry.CHROM].overlaps(start, end):
        return None, None

    cand_intervals = lookup[entry.CHROM].search(start, end)
    s_ret = min([x.data for x in cand_intervals if overlaps(start, end, x[0], x[1])])
    e_ret = max([x.data for x in cand_intervals if overlaps(start, end, x[0], x[1])])
    return s_ret, e_ret


def get_vcf_boundaries(entry):
    """
    Get the start/end of an entry and order (start < end)
    """
    start = entry.start
    if "END" in entry.INFO:
        end = entry.INFO["END"]
    else:
        end = entry.end
    if start > end:
        start, end = end, start
    return start, end

def get_vcf_entry_size(entry):
    """
    Calculate the size of the variant. Use SVLEN INFO tag if available. Otherwise inferr
    """
    if "SVLEN" in entry.INFO:
        if type(entry.INFO["SVLEN"]) is list:
            size = abs(entry.INFO["SVLEN"][0])
        else:
            size = abs(entry.INFO["SVLEN"])
    elif str(entry.ALT[0]).count("<"):
        start, end = get_vcf_boundaries(entry)
        size = end - start
    else:
        size = abs(len(entry.REF) - len(str(entry.ALT[0])))
    return size

def get_rec_ovl(astart, aend, bstart, bend):
    """
    Compute reciprocal overlap between two spans
    """
    ovl_start = max(astart, bstart)
    ovl_end = min(aend, bend)
    if ovl_start < ovl_end:  # Otherwise, they're not overlapping
        ovl_pct = float(ovl_end - ovl_start) / max(aend - astart, bend - bstart)
    else:
        ovl_pct = 0
    return ovl_pct

def get_weighted_score(sim, size, ovl):
    """
    Unite the similarity measures and make a score
    return (2*sim + 1*size + 1*ovl) / 3.0
    """
    return (2*sim + 1*size + 1*ovl) / 3.0

def edit_header(my_vcf):
    """
    Add INFO for new fields to vcf
    #Probably want to put in the PG whatever, too
    """
    # Update header
    # Edit Header
    my_vcf.infos['TruScore'] = vcf.parser._Info(
        id="TruScore", num=1, type='Float',
        desc="Truvari score for similarity of match",
        source=None, version=None)
    my_vcf.infos['PctSeqSimilarity'] = vcf.parser._Info(
        id="PctSeqSimilarity", num=1, type='Float',
        desc="Pct sequence similarity between this variant and it's closest match",
        source=None, version=None)
    my_vcf.infos['PctSizeSimilarity'] = vcf.parser._Info(
        id="PctSizeSimilarity", num=1, type='Float',
        desc="Pct size similarity between this variant and it's closest match",
        source=None, version=None)
    my_vcf.infos['StartDistance'] = vcf.parser._Info(
        id="StartDistance", num=1, type='Integer',
        desc="Distance of this call's start from comparison call's start",
        source=None, version=None)
    my_vcf.infos['EndDistance'] = vcf.parser._Info(
        id="EndDistance", num=1, type='Integer',
        desc="Distance of this call's start from comparison call's start",
        source=None, version=None)
    my_vcf.infos['PctRecOverlap'] = vcf.parser._Info(
        id="PctRecOverlap", num=1, type='Float',
        desc="Percent reciprocal overlap percent of the two calls' coordinates",
        source=None, version=None)
    my_vcf.infos['SizeDiff'] = vcf.parser._Info(
        id="SizeDiff", num=1, type='Float',
        desc="Difference in size(basecall) and size(evalcall)",
        source=None, version=None)
    my_vcf.infos['NumNeighbors'] = vcf.parser._Info(
        id="NumNeighbors", num=1, type='Integer',
        desc="Number of calls in B that were in the neighborhood (REFDIST) of this call",
        source=None, version=None)
    my_vcf.infos['NumThresholdNeighbors'] = vcf.parser._Info(
        id="NumThresholdNeighbors", num=1, type='Integer',
        desc="Number of calls in B that are within threshold distances of this call",
        source=None, version=None)


def parse_args(args):
    """
    Pull the command line parameters
    """
    def restricted_float(x):
        x = float(x)
        if x < 0.0 or x > 1.0:
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
        return x

    parser = argparse.ArgumentParser(prog="truvari", description=USAGE,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--base", type=str, required=True,
                        help="Baseline truth-set calls")
    parser.add_argument("-c", "--calls", type=str, required=True,
                        help="Comparison set of calls")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output directory")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Verbose logging")

    thresg = parser.add_argument_group("Comparison Threshold Arguments")
    thresg.add_argument("-r", "--refdist", type=int, default=500,
                        help="Max reference location distance (%(default)s)")
    thresg.add_argument("-p", "--pctsim", type=restricted_float, default=0.85,
                        help="Min percent allele sequence similarity. Set to 0 to ignore. (%(default)s)")
    thresg.add_argument("-P", "--pctsize", type=restricted_float, default=0.70,
                        help="Min pct allele size similarity (minvarsize/maxvarsize) (%(default)s)")
    #This is dumb, isn't it
    thresg.add_argument("-t", "--typeignore", action="store_true", default=False,
                        help="Compare variants of different types (%(default)s)")
    thresg.add_argument("-s", "--sizemin", type=int, default=50,
                        help="Minimum variant size to consider for comparison (%(default)s)")
    thresg.add_argument("-S", "--sizefilt", type=int, default=30,
                        help="Minimum variant size to load into IntervalTree (%(default)s)")
    thresg.add_argument("--sizemax", type=int, default=50000,
                        help="Maximum variant size to consider for comparison (%(default)s)")
    thresg.add_argument("--use-swalign", action="store_true", default=False,
                        help=("Use SmithWaterman to compute sequence similarity "
                              "instead of Levenshtein ratio. WARNING - much slower. (%(default)s)"))
    #thresg.add_argument("--use-levratio", action="store_true", default=False,
                        #help=("Use Levenshtein ratio to compute estimate similarity "
                        #"instead of editdistance. (%(default)s)"))
    filteg = parser.add_argument_group("Filtering Arguments")
    filteg.add_argument("--passonly", action="store_true", default=False,
                        help="Only consider calls with FILTER == PASS")
    genoty = parser.add_argument_group("Genotype Comparison Arguments")
    genoty.add_argument("--gtcomp", action="store_true", default=False,
                        help="Compare GenoTypes of matching calls")
    genoty.add_argument("--bSample", type=str, default=None,
                        help="Baseline calls sample to use (first)")
    genoty.add_argument("--cSample", type=str, default=None,
                        help="Comparison calls sample to use (first)")
    genoty.add_argument("--no-ref", action="store_true", default=False,
                        help="Don't include 0/0 or ./. GT calls (%(default)s)")
    args = parser.parse_args(args)
    return args


def run(cmdargs):
    args = parse_args(cmdargs)
    if os.path.isdir(args.output):
        print("Error! Output Directory %s already exists" % args.output)
        exit(1)

    os.mkdir(args.output)

    setup_logging(args.debug, LogFileStderr(os.path.join(args.output, "log.txt")))
    logging.info("Params:\n%s", json.dumps(vars(args), indent=4))

    ref_gts = ["0/0", "0|0", "./.", ".|."]

    vcfA = vcf.Reader(open(args.base, 'r'))
    if args.bSample is not None:
        sampleA = args.bSample
        if sampleA not in vcfA.samples:
            logging.error("Sample %s not found in vcf (%s)", sampleA, vcfA.samples)
            exit(1)
    else:
        sampleA = vcfA.samples[0]

    vcfB = vcf.Reader(open(args.calls, 'r'))
    edit_header(vcfB)
    if args.cSample is not None:
        sampleB = args.cSample
        if sampleB not in vcfB.samples:
            logging.error("Sample %s not found in vcf (%s)", sampleB, vcfB.samples)
            exit(1)
    else:
        sampleB = vcfB.samples[0]

    logging.info("Creating call interval tree for overlap search")
    span_lookup, num_entries_b, cmp_entries= make_interval_tree(vcfB, args.sizefilt, args.passonly)

    logging.info("%d call variants", num_entries_b)
    logging.info("%d call variants within size range (%d, %d)", cmp_entries, args.sizefilt, args.sizemax)

    num_entries = 0
    for entry in vcfA:
        num_entries += 1
    logging.info("%s base variants", num_entries)
    # Reset
    vcfA = vcf.Reader(open(args.base, 'r'))
    edit_header(vcfA)

    # Make a progress bar
    bar = progressbar.ProgressBar(redirect_stdout=True, max_value=num_entries, widgets=[
        ' [', progressbar.Timer(), ' ', progressbar.Counter(), '/', str(num_entries), '] ',
        progressbar.Bar(),
        ' (', progressbar.ETA(), ') ',
    ])

    # Setup outputs
    tpb_out = vcf.Writer(open(os.path.join(args.output, "tp-base.vcf"), 'w'), vcfA)
    b_filt = vcf.Writer(open(os.path.join(args.output, "base-filter.vcf"), 'w'), vcfA)
    tpc_out = vcf.Writer(open(os.path.join(args.output, "tp-call.vcf"), 'w'), vcfB)
    c_filt = vcf.Writer(open(os.path.join(args.output, "call-filter.vcf"), 'w'), vcfB)
    fn_out = vcf.Writer(open(os.path.join(args.output, "fn.vcf"), 'w'), vcfA)
    fp_out = vcf.Writer(open(os.path.join(args.output, "fp.vcf"), 'w'), vcfB)

    stats_box = {"TP-base": 0,
                 "TP-call": 0,
                 "FP": 0,
                 "FN": 0,
                 "base cnt": 0,
                 "call cnt": 0,
                 "base size filtered": 0,
                 "call size filtered": 0,
                 "base gt filtered": 0,
                 "call gt filtered": 0}

    # Only match calls once
    already_considered = defaultdict(bool)
    # for variant A - do var_match in B
    logging.info("Matching base to calls")
    for pbarcnt, entryA in enumerate(vcfA):
        bar.update(pbarcnt)
        sz = get_vcf_entry_size(entryA)
        
        if sz < args.sizemin or sz > args.sizemax:
            stats_box["base size filtered"] += 1
            b_filt.write_record(entryA)
            continue
        if args.no_ref and entryA.genotype(sampleA)["GT"] in ref_gts:
            stats_box["base gt filtered"] += 1
            b_filt.write_record(entryA)
            continue

        if args.passonly and len(entryA.FILTER):
            continue
        stats_box["base cnt"] += 1

        fetch_start, fetch_end = fetch_coords(span_lookup, entryA, args.refdist)
        if fetch_start is None and fetch_end is None:
            # No overlaps, don't even bother checking
            entryA.INFO["NumNeighbors"] = 0
            entryA.INFO["NumThresholdNeighbors"] = 0
            stats_box["FN"] += 1
            fn_out.write_record(entryA)
            continue

        # IntervalTree can give boundaries past REFDIST in the case of Inversions where start>end
        # We still need to fetch on the expanded boundaries so we can test them, but
        #   we need to filter calls that otherwise shouldn't be considered
        #  see the bstart/bend below
        astart, aend = get_vcf_boundaries(entryA)
        
        thresh_neighbors = []
        num_neighbors = 0

        # methodize so we can parallelize?
        # - i'm woried the vcfRecord isn't pickleable so we'd need to collect them and then return
        # more importantly, we'd also have race conditions for when it's been used...
        # even if already_considered was atomic, you may get diff answers
        # +- 1 just to be safe because why not...
        for entryB in vcfB.fetch(entryA.CHROM, max(0, fetch_start - 1), fetch_end + 1):
            # There is a race condition here that could potentially mismatch things
            # If base1 passes matching call1 and then base2 passes matching call1
            # better, it can't use it and we mismatch
            if already_considered[vcf_to_key(entryB)]:
                continue
            nsize = get_vcf_entry_size(entryB)
            if get_vcf_entry_size(entryB) < args.sizefilt:
                continue
            
            # Double ensure OVERLAP - there's a weird edge case where fetch with
            # the interval tree can return non-overlaps
            bstart, bend = get_vcf_boundaries(entryB)
            if not overlaps(astart - args.refdist, aend + args.refdist, bstart, bend):
                continue

            if args.no_ref and entryB.genotype(sampleB)["GT"] in ref_gts:
                continue
            if args.gtcomp and entryA.genotype(sampleA)["GT"] != entryB.genotype(sampleB)["GT"]:
                continue

            num_neighbors += 1

            # Size matching/neighbor sorting here, also?
            # if the sequence similarity is the same, take the one closer in size?
            size_similarity = var_sizesim(entryA, entryB, args.typeignore)
            #logging.warning("You left in test code")
            if size_similarity < args.pctsize:
                continue
            if args.pctsim > 0:
                if args.use_swalign:
                    seq_similarity = var_pctsim_sw(entryA, entryB, args.typeignore)
                else:
                    seq_similarity = var_pctsim_lev(entryA, entryB, args.typeignore)
                #use to support this manual edit distance, but the above is better
                    #seq_similarity = var_pctsim_ed(entryA, entryB, args.typeignore)
            else:
                seq_similarity = 0
            
            ovl_pct = get_rec_ovl(astart, aend, bstart, bend)

            start_distance = astart - bstart
            end_distance = aend - bend
            logging.debug(str(entryA))
            logging.debug(str(entryB))
            logging.debug("%d %d %d", aend, bend, end_distance)
            if seq_similarity >= args.pctsim and size_similarity >= args.pctsize:
                score = get_weighted_score(seq_similarity, size_similarity, ovl_pct)
                thresh_neighbors.append((score, seq_similarity, size_similarity, ovl_pct, start_distance, end_distance, entryB))

        # reporting the best match
        entryA.INFO["NumNeighbors"] = num_neighbors
        entryA.INFO["NumThresholdNeighbors"] = len(thresh_neighbors)
        if len(thresh_neighbors) > 0:
            thresh_neighbors.sort(reverse=True)
            stats_box["TP-base"] += 1
            stats_box["TP-call"] += 1
            match_score, match_pctsim, match_pctsize, matching_ovlpct, matching_stdist, \
                matching_endist, matching_entry = thresh_neighbors[0]

            matching_entry.INFO["TruScore"] = match_score
            matching_entry.INFO["NumNeighbors"] = num_neighbors
            matching_entry.INFO["NumThresholdNeighbors"] = len(thresh_neighbors)

            # Don't use it twice
            already_considered[vcf_to_key(matching_entry)] = True

            # methodize this
            entryA.INFO["PctSeqSimilarity"] = match_pctsim
            matching_entry.INFO["PctSeqSimilarity"] = match_pctsim
            entryA.INFO["PctSizeSimilarity"] = match_pctsize
            matching_entry.INFO["PctSizeSimilarity"] = match_pctsize

            entryA.INFO["StartDistance"] = matching_stdist
            entryA.INFO["EndDistance"] = matching_endist

            matching_entry.INFO["StartDistance"] = matching_stdist
            matching_entry.INFO["EndDistance"] = matching_endist
            
            entryA.INFO["PctRecOverlap"] = matching_ovlpct
            matching_entry.INFO["PctRecOverlap"] = matching_ovlpct

            size_diff = get_vcf_entry_size(entryA) - get_vcf_entry_size(matching_entry)
            entryA.INFO["SizeDiff"] = size_diff
            matching_entry.INFO["SizeDiff"] = size_diff

            tpb_out.write_record(entryA)
            tpc_out.write_record(matching_entry)
        else:
            stats_box["FN"] += 1
            fn_out.write_record(entryA)
    bar.finish()
    do_stats_math = True
    if stats_box["TP-call"] == 0 and stats_box["FN"] == 0:
        logging.warning("No TP or FN calls in base!")
        do_stats_math = False
    else:
        logging.info("Results peek: %d TP %d FN %.2f%% Recall", stats_box["TP-call"], stats_box["FN"],
                100*(float(stats_box["TP-call"]) / (stats_box["TP-call"] + stats_box["FN"])))

    logging.info("Parsing FPs from calls")
    bar = progressbar.ProgressBar(redirect_stdout=True, max_value=num_entries_b, widgets=[
        ' [', progressbar.Timer(), ' ', progressbar.Counter(), '/', str(num_entries_b), '] ',
        progressbar.Bar(),
        ' (', progressbar.ETA(), ') ',
    ])
    vcfB = vcf.Reader(open(args.calls, 'r'))
    edit_header(vcfB)
    cnt = 0
    for entry in vcfB:
        if args.passonly and len(entry.FILTER):
            continue
        bar.update(cnt)
        cnt += 1
        stats_box["call cnt"] += 1
        if already_considered[vcf_to_key(entry)]:
            continue
        size = get_vcf_entry_size(entry)
        if size < args.sizemin:
            c_filt.write_record(entry)
            stats_box["call size filtered"] += 1
        elif args.no_ref and entry.genotype(sampleB)["GT"] in ref_gts:
            stats_box["call gt filtered"] += 1
        else:
            fp_out.write_record(entry)
            stats_box["FP"] += 1
    bar.finish()
    if do_stats_math:
        # precision
        stats_box["precision"] = float(stats_box["TP-call"]) / (stats_box["TP-call"] + stats_box["FP"])
        # recall
        stats_box["recall"] = float(stats_box["TP-call"]) / (stats_box["TP-call"] + stats_box["FN"])
    else:
        stats_box["precision"] = 0
        stats_box["recall"] = 0
    # f-measure
    neum = stats_box["recall"] * stats_box["precision"]
    denom = stats_box["recall"] + stats_box["precision"]
    if denom != 0:
        stats_box["f1"] = 2 * (neum / denom)
    else:
        stats_box["f1"] = "NaN"

    with open(os.path.join(args.output, "summary.txt"), 'w') as fout:
        fout.write(json.dumps(stats_box, indent=4))
        logging.info("Stats: %s", json.dumps(stats_box, indent=4))
    logging.info("Finished")

if __name__ == '__main__':
    run(sys.argv[1:])
