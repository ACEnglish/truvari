"""
Structural variant caller comparison tool
Given a benchmark and callset, calculate the recall/precision/f-measure
"""
import os
import sys
import copy
import json
import types
import logging
import argparse
import itertools

from functools import total_ordering
from collections import defaultdict, OrderedDict, Counter

import pysam
import numpy as np

import truvari
from truvari.giab_report import make_giabreport

shared_space = types.SimpleNamespace()


@total_ordering
class MatchResult():  # pylint: disable=too-many-instance-attributes
    """
    A base/comp match holder
    """
    __slots__ = ["base", "comp", "state", "seqsim", "sizesim", "ovlpct", "sizediff",
                 "st_dist", "ed_dist", "gt_match", "multi", "score"]

    def __init__(self):
        self.base = None
        self.comp = None
        self.seqsim = None
        self.sizesim = None
        self.ovlpct = None
        self.sizediff = None
        self.st_dist = None
        self.ed_dist = None
        self.gt_match = None
        self.state = False
        self.score = None
        self.multi = None

    def calc_score(self):
        """
        Set self.score to truvari.weighted_score
        """
        if None not in [self.seqsim, self.sizesim, self.ovlpct]:
            self.score = truvari.weighted_score(
                self.seqsim, self.sizesim, self.ovlpct)

    def __lt__(self, other):
        if self.state == other.state:
            return self.score < other.score
        return self.state < other.state

    def __eq__(self, other):
        return self.state == other.state and self.score == other.score

    def __repr__(self):
        return f'<{self.score} {self.state} {self.base.chrom}:{self.base.pos}->{self.comp.chrom}:{self.comp.pos}>'

    def __str__(self):
        return f'{self.state} {self.score} ->\n {self.base} {self.comp}'


class Matcher():
    """
    Holds matching parameters. Allows calls to be checked for filtering and matches to be made

    >>> import pysam
    >>> import truvari
    >>> mat = truvari.Matcher()
    >>> mat.params.pctsim = 0
    >>> v = pysam.VariantFile('repo_utils/test_files/input1.vcf.gz')
    >>> one = next(v); two = next(v)
    >>> mat.build_match(one, two)
    <0 False chr20:66235->chr20:68303>
    """

    def __init__(self, params=None, args=None):
        if args is not None:
            params = self.make_match_params_from_args(args)
        elif params is None:
            params = self.make_match_params()

        self.params = params

    @staticmethod
    def make_match_params():
        """
        Makes a simple namespace of matching parameters
        """
        ret = types.SimpleNamespace()
        ret.reference = None
        ret.refdist = 500
        ret.pctsim = 0.70
        ret.buffer = 0.10
        ret.pctsize = 0.70
        ret.pctovl = 0.0
        ret.typeignore = False
        ret.use_lev = False
        ret.chunksize = 1000
        ret.gtcomp = False
        ret.bSample = None
        ret.cSample = None
        # filtering properties
        ret.sizemin = 50
        ret.sizefilt = 30
        ret.sizemax = 50000
        ret.passonly = False
        ret.no_ref = False
        ret.includebed = None
        ret.multimatch = False
        return ret

    @staticmethod
    def make_match_params_from_args(args):
        """
        Makes a simple namespace of matching parameters
        """
        ret = types.SimpleNamespace()
        ret.reference = pysam.FastaFile(
            args.reference) if args.reference else None
        ret.refdist = args.refdist
        ret.pctsim = args.pctsim
        ret.buffer = args.buffer
        ret.pctsize = args.pctsize
        ret.pctovl = args.pctovl
        ret.typeignore = args.typeignore
        ret.use_lev = args.use_lev
        ret.chunksize = args.chunksize
        ret.gtcomp = args.gtcomp
        ret.bSample = args.bSample
        ret.cSample = args.cSample
        # filtering properties
        ret.sizemin = args.sizemin
        ret.sizefilt = args.sizefilt
        ret.sizemax = args.sizemax
        ret.passonly = args.passonly
        ret.no_ref = args.no_ref
        ret.multimatch = args.multimatch
        return ret

    def filter_call(self, entry, base=False):
        """
        Returns True if the call should be filtered
        Base has different filtering requirements, so let the method know
        """
        # See if it hits the regions, first
        if self.params.includebed and not self.params.includebed.include(entry):
            return True

        size = truvari.entry_size(entry)
        args = self.params
        if base and size < args.sizemin or size > args.sizemax:
            return True

        if not base and size < args.sizefilt or size > args.sizemax:
            return True

        samp = args.bSample if base else args.cSample
        prefix = 'b' if base else 'c'
        if args.no_ref in ["a", prefix] and not truvari.entry_is_present(entry, samp):
            return True

        if args.passonly and truvari.filter_value(entry):
            return True

        return False

    def build_match(self, base, comp):
        """
        Build a MatchResult
        """
        ret = MatchResult()
        ret.base = base
        ret.comp = comp
        ret.state = True

        if not self.params.typeignore and not truvari.entry_same_variant_type(base, comp):
            logging.debug("%s and %s are not the same SVTYPE",
                          str(base), str(comp))
            ret.state = False

        bstart, bend = truvari.entry_boundaries(base)
        cstart, cend = truvari.entry_boundaries(comp)
        if not truvari.overlaps(bstart - self.params.refdist, bend + self.params.refdist, cstart, cend):
            logging.debug("%s and %s are not within REFDIST",
                          str(base), str(comp))
            ret.state = False

        ret.sizesim, ret.sizediff = truvari.entry_size_similarity(base, comp)
        if ret.sizesim < self.params.pctsize:
            logging.debug("%s and %s size similarity is too low (%.3f)",
                          str(base), str(comp), ret.sizesim)
            ret.state = False

        ret.gt_match = truvari.entry_gt_comp(
            base, comp, self.params.bSample, self.params.cSample)
        if self.params.gtcomp and not ret.gt_match:
            logging.debug("%s and %s are not the same genotype",
                          str(base), str(comp))
            ret.state = False

        ret.ovlpct = truvari.entry_reciprocal_overlap(base, comp)
        if truvari.entry_variant_type(base) == "DEL" and ret.ovlpct < self.params.pctovl:
            logging.debug("%s and %s overlap percent is too low (%.3f)",
                          str(base), str(comp), ret.ovlpct)
            ret.state = False

        if self.params.pctsim > 0:
            ret.seqsim = truvari.entry_pctsim(base, comp, self.params.reference,
                                              self.params.buffer, self.params.use_lev)
            if ret.seqsim < self.params.pctsim:
                logging.debug("%s and %s sequence similarity is too low (%.3ff)",
                              str(base), str(comp), ret.seqsim)
                ret.state = False
        else:
            ret.seqsim = 0

        ret.st_dist, ret.ed_dist = truvari.entry_distance(base, comp)
        ret.calc_score()

        return ret


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
        logging.info("Results peek: %d TP-base %d FN %.2f%% Recall", self["TP-base"], self["FN"],
                     100 * (float(self["TP-base"]) / (self["TP-base"] + self["FN"])))

        # Final calculations
        if do_stats_math:
            self["precision"] = float(
                self["TP-call"]) / (self["TP-call"] + self["FP"])
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
        else:
            self["f1"] = "NaN"


############################
# Parsing and set building #
############################
def file_zipper(*start_files):
    """
    Zip files to yield the entries in order.
    Files must be sorted
    Each parameter is a tuple of ('key', 'fn')
    where key is the identifier (so we know which file the yielded entry came from)
    and fn is a VariantFile
    yields key, pysam.VariantRecord
    """
    next_markers = []
    files = []
    names = []
    file_counts = Counter()
    for name, i in start_files:
        try:
            next_markers.append(next(i))
            files.append(i)
            names.append(name)
        except StopIteration:
            # For when there are no variants in the file
            pass

    while next_markers:
        sidx = 0 # assume the first is the least
        for idx, i in enumerate(next_markers):
            if i.chrom < next_markers[sidx].chrom:
                sidx = idx
            elif i.chrom == next_markers[sidx].chrom and i.start < next_markers[sidx].start:
                sidx = idx
        entry = next_markers[sidx]
        key = names[sidx]
        file_counts[key] += 1
        try:
            next_markers[sidx] = next(files[sidx])
        except StopIteration:
            # This file is done
            files.pop(sidx)
            names.pop(sidx)
            next_markers.pop(sidx)
        yield key, entry
    logging.info(f"Zipped {sum(file_counts.values())} variants. {file_counts}")


def chunker(matcher, *files):
    """
    Given a Matcher and multiple files,
    """
    cur_chrom = None
    #min_start = None
    max_end = None
    call_counts = Counter()
    chunk_count = 0
    cur_chunk = defaultdict(list)
    for key, entry in file_zipper(*files):
        new_chunk = max_end and max_end + matcher.params.chunksize < entry.start
        new_chrom = cur_chrom and entry.chrom != cur_chrom
        if new_chunk or new_chrom:
            chunk_count += 1
            yield matcher, cur_chunk
            cur_chunk = defaultdict(list)
            cur_chrom = None
            #min_start = None
            max_end = None  # Gotta encounter a new SV to get a new prev_pos

        if not matcher.filter_call(entry, key == 'base'):
            if not cur_chrom:  # First call in chunk
                cur_chrom = entry.chrom
                #min_start = entry.start
            max_end = entry.stop
            logging.debug(f"Adding to {key} -> {entry}")
            cur_chunk[key].append(entry)
            call_counts[key] += 1
    chunk_count += 1
    logging.info(f"{chunk_count} chunks of {sum(call_counts.values())} variants. {call_counts}")
    return matcher, cur_chunk


def compare_chunk(chunk):
    """
    Given a filtered chunk, (from chunker) compare all of the calls
    """
    # If I'm giving up on multiprocessing (for now) I can move back to making the chunker do the holding
    # Instead of fetching. Then I can get rid of base_vcf and comp_vcf in shared_space
    matcher, chunk_dict = chunk
    logging.debug(f"Comparing chunk {chunk_dict}")
    # Build all v all matrix
    match_matrix = []
    fns = []
    for b in chunk_dict['base']:
        base_matches = []
        for c in chunk_dict['comp']:
            mat = matcher.build_match(b, c)
            logging.debug(f"Made mat -> {mat}")
            base_matches.append(mat)
        if not base_matches:  # All FNs
            ret = MatchResult()
            ret.base = b
            logging.debug(f"All FN chunk -> {ret}")
            fns.append(ret)
        else:
            match_matrix.append(base_matches)

    # need to iterate the comp again and make them all FPs
    if not match_matrix:
        fps = []
        for c in chunk_dict['comp']:
            ret = MatchResult()
            ret.comp = c
            fps.append(ret)
            logging.debug("All FP chunk -> {ret}")
        return fps + fns

    # Sort and annotate matches
    match_matrix = np.array(match_matrix)
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
        if base_cnt == 0 and not comp_is_used:  # Only pull the comp
            to_process = copy.copy(match)
            to_process.base = None
            to_process.state = False
            to_process.multi = True
            comp_cnt -= 1
            used_comp.add(str(to_process.comp))
            ret.append(to_process)
        elif comp_cnt == 0 and not base_is_used:  # Only pull the base
            to_process = copy.copy(match)
            to_process.comp = None
            to_process.state = False
            to_process.multi = True
            base_cnt -= 1
            used_base.add(str(to_process.base))
            ret.append(to_process)
        # A match
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
    #Probably want to put in the PG whatever, too
    """
    # Update header
    # Edit Header
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
                     'Description="Distance of this call\'s start from comparison call\'s start">'))
    header.add_line(('##INFO=<ID=EndDistance,Number=1,Type=Integer,'
                     'Description="Distance of this call\'s start from comparison call\'s start">'))
    header.add_line(('##INFO=<ID=SizeDiff,Number=1,Type=Float,'
                     'Description="Difference in size(basecall) and size(evalcall)">'))
    header.add_line(('##INFO=<ID=GTMatch,Number=0,Type=Flag,'
                     'Description="Base/Comparison Genotypes match">'))
    header.add_line(('##INFO=<ID=MatchId,Number=1,Type=Integer,'
                     'Description="Truvari uid to help tie tp-base.vcf and tp-call.vcf entries together">'))
    header.add_line(('##INFO=<ID=Multi,Number=0,Type=Flag,'
                     'Description="Call is false due to non-multimatching">'))
    return header


def annotate_entry(entry, match, header):
    """
    Make a new entry with all the information
    """
    entry = truvari.copy_entry(entry, header)
    entry.info["PctSeqSimilarity"] = match.seqsim
    entry.info["PctSizeSimilarity"] = match.sizesim
    entry.info["PctRecOverlap"] = match.ovlpct
    entry.info["SizeDiff"] = match.sizediff
    entry.info["StartDistance"] = match.st_dist
    entry.info["EndDistance"] = match.ed_dist
    entry.info["GTMatch"] = match.gt_match
    entry.info["TruScore"] = match.score
    entry.info["Multi"] = match.multi
    return entry


def output_writer(call):
    """
    Given a call or match result or something
    write it to the appropriate file
    This will need to figure out where its going
    Can also call annotate here, I reckon
    And I'll deal with the stats box
    """
    outs = shared_space.outputs
    box = shared_space.outputs["stats_box"]
    if call.base:
        box["base cnt"] += 1
        entry = annotate_entry(call.base, call, outs['n_base_header'])
        if call.state:
            box["TP-base"] += 1
            outs["tpb_out"].write(entry)
            if call.gt_match:
                box["TP-base_TP-gt"] += 1
            else:
                box["TP-base_FP-gt"] += 1
        else:
            box["FN"] += 1
            outs["fn_out"].write(entry)
    if call.comp:
        entry = annotate_entry(call.comp, call, outs['n_comp_header'])
        if call.state:
            box["call cnt"] += 1
            box["TP-call"] += 1
            outs["tpc_out"].write(entry)
            if call.gt_match:
                box["TP-call_TP-gt"] += 1
            else:
                box["TP-call_FP-gt"] += 1
        elif truvari.entry_size(entry) >= shared_space.sizemin:
            # The if is because we don't count FPs between sizefilt-sizemin
            box["call cnt"] += 1
            box["FP"] += 1
            outs["fp_out"].write(entry)


##########################
# Parameters and outputs #
##########################
def parse_args(args):
    """
    Pull the command line parameters
    """
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
    thresg.add_argument("-r", "--refdist", type=int, default=500,
                        help="Max reference location distance (%(default)s)")
    thresg.add_argument("-p", "--pctsim", type=truvari.restricted_float, default=0.70,
                        help="Min percent allele sequence similarity. Set to 0 to ignore. (%(default)s)")
    thresg.add_argument("-B", "--buffer", type=truvari.restricted_float, default=0.10,
                        help="Percent of the reference span to buffer the haplotype sequence created")
    thresg.add_argument("-P", "--pctsize", type=truvari.restricted_float, default=0.70,
                        help="Min pct allele size similarity (minvarsize/maxvarsize) (%(default)s)")
    thresg.add_argument("-O", "--pctovl", type=truvari.restricted_float, default=0.0,
                        help="Min pct reciprocal overlap (%(default)s) for DEL events")
    thresg.add_argument("-t", "--typeignore", action="store_true", default=False,
                        help="Variant types don't need to match to compare (%(default)s)")
    thresg.add_argument("--use-lev", action="store_true",
                        help="Use the Levenshtein distance ratio instead of edlib editDistance ratio (%(default)s)")
    thresg.add_argument("-C", "--chunksize", type=int, default=None,
                        help="Max reference distance to compare calls (2 x refdist)")

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
    if args.chunksize is None:
        args.chunksize = 2 * args.refdist
    if args.chunksize < args.refdist:
        parser.error("--chunksize must be >= --refdist")
    return args


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
        logging.error(
            "Comparison vcf %s does not end with .gz. Must be bgzip'd", args.comp)
    if not os.path.exists(args.comp + '.tbi'):
        check_fail = True
        logging.error(
            "Comparison vcf index %s.tbi does not exist. Must be indexed", args.comp)
    if not args.base.endswith(".gz"):
        check_fail = True
        logging.error(
            "Base vcf %s does not end with .gz. Must be bgzip'd", args.base)
    if not os.path.exists(args.base + '.tbi'):
        check_fail = True
        logging.error(
            "Base vcf index %s.tbi does not exist. Must be indexed", args.base)

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
    Places the data into the shared space
    """
    os.mkdir(args.output)
    truvari.setup_logging(args.debug, truvari.LogFileStderr(
        os.path.join(args.output, "log.txt")))
    logging.info("Params:\n%s", json.dumps(vars(args), indent=4))
    logging.info(f"Truvari version: {truvari.__version__}")
    shared_space.outputs = {}

    shared_space.outputs["vcf_base"] = pysam.VariantFile(args.base)
    shared_space.outputs["n_base_header"] = edit_header(
        shared_space.outputs["vcf_base"])
    shared_space.outputs["sampleBase"] = args.bSample if args.bSample else shared_space.outputs["vcf_base"].header.samples[0]

    shared_space.outputs["vcf_comp"] = pysam.VariantFile(args.comp)
    shared_space.outputs["n_comp_header"] = edit_header(
        shared_space.outputs["vcf_comp"])
    shared_space.outputs["sampleComp"] = args.cSample if args.cSample else shared_space.outputs["vcf_comp"].header.samples[0]

    # Setup outputs
    shared_space.outputs["tpb_out"] = pysam.VariantFile(os.path.join(
        args.output, "tp-base.vcf"), 'w', header=shared_space.outputs["n_base_header"])
    shared_space.outputs["tpc_out"] = pysam.VariantFile(os.path.join(
        args.output, "tp-call.vcf"), 'w', header=shared_space.outputs["n_comp_header"])

    shared_space.outputs["fn_out"] = pysam.VariantFile(os.path.join(
        args.output, "fn.vcf"), 'w', header=shared_space.outputs["n_base_header"])
    shared_space.outputs["fp_out"] = pysam.VariantFile(os.path.join(
        args.output, "fp.vcf"), 'w', header=shared_space.outputs["n_comp_header"])

    shared_space.outputs["stats_box"] = StatsBox()


def close_outputs():
    """
    Close all the files
    """
    shared_space.outputs["tpb_out"].close()
    shared_space.outputs["tpc_out"].close()
    shared_space.outputs["fn_out"].close()
    shared_space.outputs["fp_out"].close()


def bench_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)

    # Need to build the MatchParams
    # and I guess I can still use shared_space?? For outputs, mainly
    # So let's contain the outputs to just a section that's specific
    # to bench.py only
    if check_params(args) or check_inputs(args):
        sys.stderr.write("Couldn't run Truvari. Please fix parameters\n")
        sys.exit(100)

    # build the MatchParams
    if args.reference:
        shared_space.reference = pysam.Fastafile(args.reference)

    # build_match needs to not make FPs less than sizefilt
    shared_space.sizemin = args.sizemin
    matcher = Matcher(args=args)
    setup_outputs(args)

    base = pysam.VariantFile(args.base)
    comp = pysam.VariantFile(args.comp)
    matcher.params.includebed = truvari.GenomeTree(base, comp,
                                                   args.includebed,
                                                   args.sizemax)

    chunks = chunker(matcher, ('base', base), ('comp', comp))
    for call in itertools.chain.from_iterable(map(compare_chunk, chunks)):
        output_writer(call)

    with open(os.path.join(args.output, "summary.txt"), 'w') as fout:
        box = shared_space.outputs["stats_box"]
        box.calc_performance()
        fout.write(json.dumps(box, indent=4))
        logging.info("Stats: %s", json.dumps(box, indent=4))

    close_outputs()

    if args.giabreport:
        make_giabreport(args, shared_space.outputs["stats_box"])

    logging.info("Finished bench")
