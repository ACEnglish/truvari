"""
Structural variant benchmarking tool
Given baseline and comparison sets of variants, calculate the recall/precision/f-measure
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

##############
# Parameters #
##############


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
                        help="Fasta used to call variants. Turns on reference context sequence comparison")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Verbose logging")

    thresg = parser.add_argument_group("Comparison Threshold Arguments")
    thresg.add_argument("-r", "--refdist", type=truvari.restricted_int, default=defaults.refdist,
                        help="Max reference distance (%(default)s)")
    thresg.add_argument("-p", "--pctseq", type=truvari.restricted_float, default=defaults.pctseq,
                        help="Min sequence similarity. Set to 0 to ignore (%(default)s)")
    thresg.add_argument("-P", "--pctsize", type=truvari.restricted_float, default=defaults.pctsize,
                        help="Min variant size similarity (%(default)s)")
    thresg.add_argument("-O", "--pctovl", type=truvari.restricted_float, default=defaults.pctovl,
                        help="Min reciprocal overlap (%(default)s)")
    thresg.add_argument("-t", "--typeignore", action="store_true", default=defaults.typeignore,
                        help="Don't compare variant types (%(default)s)")
    thresg.add_argument("--pick", type=str, default=defaults.pick, choices=PICKERS.keys(),
                        help="Number of matches reported per-call (%(default)s)")
    thresg.add_argument("--dup-to-ins", action="store_true",
                        help="Assume DUP svtypes are INS (%(default)s)")
    thresg.add_argument("-C", "--chunksize", type=truvari.restricted_int, default=defaults.chunksize,
                        help="Max reference distance to compare calls (%(default)s)")
    thresg.add_argument("-B", "--minhaplen", type=truvari.restricted_int, default=defaults.minhaplen,
                        help="Min haplotype sequence length to create (%(default)s)")

    genoty = parser.add_argument_group("Genotype Comparison Arguments")
    genoty.add_argument("--bSample", type=str, default=None,
                        help="Baseline calls' sample to use (first)")
    genoty.add_argument("--cSample", type=str, default=None,
                        help="Comparison calls' sample to use (first)")

    filteg = parser.add_argument_group("Filtering Arguments")
    filteg.add_argument("--passonly", action="store_true", default=defaults.passonly,
                        help="Only consider calls with FILTER == PASS")
    filteg.add_argument("-s", "--sizemin", type=truvari.restricted_int, default=defaults.sizemin,
                        help="Minimum variant size to consider from --comp (%(default)s)")
    filteg.add_argument("-S", "--sizefilt", type=truvari.restricted_int, default=None,
                        help="Minimum variant size to consider from --base (30)")
    filteg.add_argument("--sizemax", type=truvari.restricted_int, default=defaults.sizemax,
                        help="Maximum variant size to consider (%(default)s)")
    filteg.add_argument("--no-ref", default=defaults.no_ref, choices=['a', 'b', 'c'],
                        help="Exclude 0/0 or ./. GT calls from all (a), base (b), or comp (c) vcfs (%(default)s)")
    filteg.add_argument("--includebed", type=str, default=None,
                        help="Bed file of regions to analyze. Only calls within regions are counted")
    filteg.add_argument("--extend", type=truvari.restricted_int, default=0,
                        help="Distance to allow comp entries outside of includebed regions (%(default)s)")

    args = parser.parse_args(args)
    # When sizefilt is not provided and sizemin has been lowered below the default,
    # set sizefilt to sizemin. Otherwise, if sizefilt not provided, set to default
    # This just makes it easier to specify e.g. `-s 1` instead of `-s 1 -S 1`
    if args.sizefilt is None:
        if args.sizemin < defaults.sizefilt:
            args.sizefilt = args.sizemin
        else:
            args.sizefilt = defaults.sizefilt

    # Setup abspaths
    args.base = os.path.abspath(args.base)
    args.comp = os.path.abspath(args.comp)
    args.includebed = os.path.abspath(
        args.includebed) if args.includebed else args.includebed
    args.reference = os.path.abspath(
        args.reference) if args.reference else args.reference

    return args


def check_params(args):
    """
    Checks parameters as much as possible.
    All errors are written to stderr without logging since failures mean no output
    """
    check_fail = False
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
        logging.error("Comparison vcf %s does not end with .gz. Must be bgzip'd",
                      args.comp)
        check_fail = True
    if not os.path.exists(args.comp + '.tbi'):
        logging.error("Comparison vcf index %s.tbi does not exist. Must be indexed",
                      args.comp)
        check_fail = True
    if not args.base.endswith(".gz"):
        logging.error("Base vcf %s does not end with .gz. Must be bgzip'd",
                      args.base)
        check_fail = True
    if not os.path.exists(args.base + '.tbi'):
        logging.error("Base vcf index %s.tbi does not exist. Must be indexed",
                      args.base)
        check_fail = True
    if args.includebed and not os.path.exists(args.includebed):
        logging.error("Include bed %s does not exist", args.includebed)
        check_fail = True
    if args.reference and not os.path.exists(args.reference):
        logging.error("Reference %s does not exist", args.reference)
        check_fail = True
    return check_fail


def check_sample(vcf_fn, sample_id=None):
    """
    Checks that a sample is inside a vcf
    Returns True if check failed
    """
    vcf_file = pysam.VariantFile(vcf_fn)
    check_fail = False
    if sample_id is not None and sample_id not in vcf_file.header.samples:
        logging.error("Sample %s not found in vcf (%s)", sample_id, vcf_fn)
        check_fail = True
    if len(vcf_file.header.samples) == 0:
        logging.error("No SAMPLE columns found in vcf (%s)", vcf_fn)
        check_fail = True
    if sample_id is None:
        sample_id = vcf_file.header.samples[0]
    return check_fail, sample_id


def check_inputs(args):
    """
    Checks the inputs to ensure expected values are found inside of files
    Returns True if check failed
    """
    b_check, args.bSample = check_sample(args.base, args.bSample)
    c_check, args.cSample = check_sample(args.comp, args.cSample)
    return b_check or c_check


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
    header.add_line(('##INFO=<ID=GTMatch,Number=1,Type=Integer,'
                     'Description="Base/Comparison genotypes AC difference">'))
    header.add_line(('##INFO=<ID=MatchId,Number=.,Type=String,'
                     'Description="Tuple of base and comparison call ids which were matched">'))
    header.add_line(('##INFO=<ID=Multi,Number=0,Type=Flag,'
                     'Description="Call is false due to non-multimatching">'))
    return header


def annotate_entry(entry, match, header):
    """
    Make a new entry with all the information
    """
    entry.translate(header)
    entry.info["PctSeqSimilarity"] = round(
        match.seqsim, 4) if match.seqsim is not None else None
    entry.info["PctSizeSimilarity"] = round(
        match.sizesim, 4) if match.sizesim is not None else None
    entry.info["PctRecOverlap"] = round(
        match.ovlpct, 4) if match.ovlpct is not None else None
    entry.info["SizeDiff"] = match.sizediff
    entry.info["StartDistance"] = match.st_dist
    entry.info["EndDistance"] = match.ed_dist
    entry.info["GTMatch"] = match.gt_match
    entry.info["TruScore"] = int(match.score) if match.score else None
    entry.info["MatchId"] = match.matid
    entry.info["Multi"] = match.multi


#############
# Core code #
#############
class StatsBox(OrderedDict):
    """
    Make a blank stats box for counting TP/FP/FN and calculating performance
    """

    def __init__(self):
        super().__init__()
        self["TP-base"] = 0
        self["TP-comp"] = 0
        self["FP"] = 0
        self["FN"] = 0
        self["precision"] = 0
        self["recall"] = 0
        self["f1"] = 0
        self["base cnt"] = 0
        self["comp cnt"] = 0
        self["TP-comp_TP-gt"] = 0
        self["TP-comp_FP-gt"] = 0
        self["TP-base_TP-gt"] = 0
        self["TP-base_FP-gt"] = 0
        self["gt_concordance"] = 0
        self["gt_matrix"] = defaultdict(Counter)

    def calc_performance(self):
        """
        Calculate the precision/recall
        """
        if self["TP-base"] == 0 and self["FN"] == 0:
            logging.warning("No TP or FN calls in base!")
        elif self["TP-comp"] == 0 and self["FP"] == 0:
            logging.warning("No TP or FP calls in comp!")

        precision, recall, f1 = truvari.performance_metrics(
            self["TP-base"], self["TP-comp"], self["FN"], self["FP"])

        self["precision"] = precision
        self["recall"] = recall
        self["f1"] = f1
        if self["TP-comp_TP-gt"] + self["TP-comp_FP-gt"] != 0:
            self["gt_concordance"] = float(self["TP-comp_TP-gt"]) / (self["TP-comp_TP-gt"] +
                                                                     self["TP-comp_FP-gt"])

    def clean_out(self):
        """
        When reusing a StatsBox (typically inside refine), gt numbers
        are typically invalidated. This removes those numbers from self to make
        a cleaner report
        """
        del self["TP-comp_TP-gt"]
        del self["TP-comp_FP-gt"]
        del self["TP-base_TP-gt"]
        del self["TP-base_FP-gt"]
        del self["gt_concordance"]
        del self["gt_matrix"]

    def write_json(self, out_name):
        """
        Write stats as json to file
        """
        with open(out_name, 'w') as fout:
            fout.write(json.dumps(self, indent=4))


class BenchOutput():
    """
    Makes all of the output files for a Bench.run

    The variable `BenchOutput.vcf_filenames` holds a dictonary. The keys are tpb, tpc
    for true positive base/comp vcf filename and fn, fp. The variable `stats_box` holds
    a :class:`StatsBox`.
    """

    def __init__(self, bench, matcher):
        """
        initialize
        """
        self.m_bench = bench
        self.m_matcher = matcher

        os.mkdir(self.m_bench.outdir)
        param_dict = self.m_bench.param_dict()
        param_dict.update(vars(self.m_matcher.params))

        if self.m_bench.do_logging:
            truvari.setup_logging(self.m_bench.debug, truvari.LogFileStderr(
                os.path.join(self.m_bench.outdir, "log.txt")), show_version=True)
            logging.info("Params:\n%s", json.dumps(param_dict, indent=4))

        with open(os.path.join(self.m_bench.outdir, 'params.json'), 'w') as fout:
            json.dump(param_dict, fout)

        b_vcf = pysam.VariantFile(self.m_bench.base_vcf)
        c_vcf = pysam.VariantFile(self.m_bench.comp_vcf)
        self.n_headers = {'b': edit_header(b_vcf),
                          'c': edit_header(c_vcf)}

        self.vcf_filenames = {'tpb': os.path.join(self.m_bench.outdir, "tp-base.vcf"),
                              'tpc': os.path.join(self.m_bench.outdir, "tp-comp.vcf"),
                              'fn': os.path.join(self.m_bench.outdir, "fn.vcf"),
                              'fp': os.path.join(self.m_bench.outdir, "fp.vcf")}
        self.out_vcfs = {}
        for key in ['tpb', 'fn']:
            self.out_vcfs[key] = pysam.VariantFile(
                self.vcf_filenames[key], mode='w', header=self.n_headers['b'])
        for key in ['tpc', 'fp']:
            self.out_vcfs[key] = pysam.VariantFile(
                self.vcf_filenames[key], mode='w', header=self.n_headers['c'])

        self.stats_box = StatsBox()

    def write_match(self, match):
        """
        Annotate a MatchResults' entries then write to the apppropriate file
        and do the stats counting.
        Writer is responsible for handling FPs between sizefilt-sizemin
        """
        box = self.stats_box
        if match.base:
            box["base cnt"] += 1
            annotate_entry(match.base, match, self.n_headers['b'])
            if match.state:
                gtBase = str(match.base_gt)
                gtComp = str(match.comp_gt)
                box["gt_matrix"][gtBase][gtComp] += 1

                box["TP-base"] += 1
                self.out_vcfs["tpb"].write(match.base)
                if match.gt_match == 0:
                    box["TP-base_TP-gt"] += 1
                else:
                    box["TP-base_FP-gt"] += 1
            else:
                box["FN"] += 1
                self.out_vcfs["fn"].write(match.base)

        if match.comp:
            annotate_entry(match.comp, match, self.n_headers['c'])
            if match.state:
                box["comp cnt"] += 1
                box["TP-comp"] += 1
                self.out_vcfs["tpc"].write(match.comp)
                if match.gt_match == 0:
                    box["TP-comp_TP-gt"] += 1
                else:
                    box["TP-comp_FP-gt"] += 1
            elif truvari.entry_size(match.comp) >= self.m_matcher.params.sizemin:
                # The if is because we don't count FPs between sizefilt-sizemin
                box["comp cnt"] += 1
                box["FP"] += 1
                self.out_vcfs["fp"].write(match.comp)



    def close_outputs(self):
        """
        Close all the files
        """
        for i in self.out_vcfs.values():
            i.close()

        for i in self.vcf_filenames.values():
            truvari.compress_index_vcf(i)

        self.stats_box.calc_performance()
        self.stats_box.write_json(os.path.join(
            self.m_bench.outdir, "summary.json"))


class Bench():
    """
    Object to perform operations of truvari bench

    You can build a Bench worker with default matching parameters via:

    .. code-block:: python

       m_bench = truvari.Bench()

    If you'd like to customize the parameters, build and edit a :class:`Matcher` and pass it to the Bench init

    .. code-block:: python

       matcher = truvari.Matcher()
       matcher.params.pctseq = 0.50
       m_bench = truvari.Bench(matcher)

    To run on a chunk of :class:`pysam.VariantRecord` already loaded, pass them in as lists to:

    .. code-block:: python

       matches = m_bench.compare_calls(base_variants, comp_variants)

    For advanced parsing, one can build the match matrix for a chunk of calls with:

    .. code-block:: python

       match_matrix = m_bench.build_matrix(base_variants, comp_variants)

    This can then be used for things such as creating custom match pickers or adding new matching checks.

    If you want to run on existing files, simply supply the arguments to init

    .. code-block:: python

        m_bench = truvari.Bench(matcher, base_vcf, comp_vcf, outdir)
        output = m_bench.run()

    Note that running on files must write to an output directory and is the only way to use things like 'includebed'.
    However, the returned `BenchOutput` has attributes pointing to all the results.
    """

    def __init__(self, matcher=None, base_vcf=None, comp_vcf=None, outdir=None,
                 includebed=None, extend=0, debug=False, do_logging=False):
        """
        Initilize
        """
        self.matcher = matcher if matcher is not None else truvari.Matcher()
        self.base_vcf = base_vcf
        self.comp_vcf = comp_vcf
        self.outdir = outdir
        self.includebed = includebed
        self.extend = extend
        self.debug = debug
        self.do_logging = do_logging
        self.refine_candidates = []

    def param_dict(self):
        """
        Returns the parameters as a dict
        """
        return {'base': self.base_vcf,
                "comp": self.comp_vcf,
                "output": self.outdir,
                "includebed": self.includebed,
                "extend": self.extend,
                "debug": self.debug}

    def run(self):
        """
        Runs bench and returns the resulting :class:`truvari.BenchOutput`
        """
        if self.base_vcf is None or self.comp_vcf is None or self.outdir is None:
            raise RuntimeError(
                "Cannot call Bench.run without base/comp vcf filenames and outdir")

        output = BenchOutput(self, self.matcher)

        base = pysam.VariantFile(self.base_vcf)
        comp = pysam.VariantFile(self.comp_vcf)

        regions = truvari.RegionVCFIterator(base, comp,
                                            self.includebed,
                                            self.matcher.params.sizemax)
        regions.merge_overlaps()
        regions_extended = regions.extend(
            self.extend) if self.extend else regions

        base_i = regions.iterate(base)
        comp_i = regions_extended.iterate(comp)

        chunks = truvari.chunker(
            self.matcher, ('base', base_i), ('comp', comp_i))
        for match in itertools.chain.from_iterable(map(self.compare_chunk, chunks)):
            # setting non-matched comp variants that are not fully contained in the original regions to None
            # These don't count as FP or TP and don't appear in the output vcf files
            if self.extend and match.comp is not None and not match.state and not regions.include(match.comp):
                match.comp = None
            output.write_match(match)

        with open(os.path.join(self.outdir, 'candidate.refine.bed'), 'w') as fout:
            fout.write("\n".join(self.refine_candidates))

        output.close_outputs()
        return output

    def compare_chunk(self, chunk):
        """
        Given a filtered chunk, (from chunker) compare all of the calls
        """
        chunk_dict, chunk_id = chunk
        logging.debug("Comparing chunk %s", chunk_id)
        result = self.compare_calls(
            chunk_dict["base"], chunk_dict["comp"], chunk_id)
        self.check_refine_candidate(result)
        return result

    def compare_calls(self, base_variants, comp_variants, chunk_id=0):
        """
        Builds MatchResults, returns them as a numpy matrix if there's at least one base and one comp variant.
        Otherwise, returns a list of the variants placed in MatchResults
        """
        # All FPs
        if len(base_variants) == 0:
            fps = []
            for cid, c in enumerate(comp_variants):
                ret = truvari.MatchResult()
                ret.comp = c
                ret.matid = ["", f"{chunk_id}.{cid}"]
                fps.append(ret)
                logging.debug("All FP -> %s", ret)
            return fps

        # All FNs
        if len(comp_variants) == 0:
            fns = []
            for bid, b in enumerate(base_variants):
                ret = truvari.MatchResult()
                ret.base = b
                ret.matid = [f"{chunk_id}.{bid}", ""]
                logging.debug("All FN -> %s", ret)
                fns.append(ret)
            return fns

        match_matrix = self.build_matrix(
            base_variants, comp_variants, chunk_id)
        if isinstance(match_matrix, list):
            return match_matrix
        return PICKERS[self.matcher.params.pick](match_matrix)

    def build_matrix(self, base_variants, comp_variants, chunk_id=0, skip_gt=False):
        """
        Builds MatchResults, returns them as a numpy matrix
        """
        if not base_variants or not comp_variants:
            raise RuntimeError(
                "Expected at least one base and one comp variant")
        match_matrix = []
        for bid, b in enumerate(base_variants):
            base_matches = []
            for cid, c in enumerate(comp_variants):
                mat = self.matcher.build_match(
                    b, c, [f"{chunk_id}.{bid}", f"{chunk_id}.{cid}"], skip_gt)
                logging.debug("Made mat -> %s", mat)
                base_matches.append(mat)
            match_matrix.append(base_matches)

        return np.array(match_matrix)

    def check_refine_candidate(self, result):
        """
        Adds this region as a candidate for refinement if there are unmatched variants
        """
        has_unmatched = False
        pos = []
        chrom = None
        for match in result:
            has_unmatched |= not match.state
            if match.base is not None:
                chrom = match.base.chrom
                pos.extend(truvari.entry_boundaries(match.base))
            if match.comp is not None:
                chrom = match.comp.chrom
                pos.extend(truvari.entry_boundaries(match.comp))
        if has_unmatched and pos:
            self.refine_candidates.append(f"{chrom}\t{min(*pos)}\t{max(*pos)}")


#################
# Match Pickers #
#################
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


def pick_ac_matches(match_matrix):
    """
    Given a numpy array of MatchResults
    Find upto allele count mumber of matches
    """
    ret = []
    base_cnt, comp_cnt = match_matrix.shape
    match_matrix = np.ravel(match_matrix)
    match_matrix.sort()
    used_comp = Counter()
    used_base = Counter()
    for match in match_matrix[::-1]:
        # No more matches to find
        if base_cnt == 0 and comp_cnt == 0:
            break
        b_key = str(match.base)
        c_key = str(match.comp)
        # This is a trick
        base_is_used = used_base[b_key] >= match.base_gt_count
        comp_is_used = used_comp[c_key] >= match.comp_gt_count
        # Only write the comp (FP)
        if base_cnt == 0 and not comp_is_used:
            to_process = copy.copy(match)
            to_process.base = None
            to_process.multi = to_process.state
            to_process.state = False
            comp_cnt -= 1
            if used_comp[c_key] == 0:  # Only write as F if it hasn't been a T
                ret.append(to_process)
            used_comp[c_key] = 9
        # Only write the base (FN)
        elif comp_cnt == 0 and not base_is_used:
            to_process = copy.copy(match)
            to_process.comp = None
            to_process.multi = to_process.state
            to_process.state = False
            base_cnt -= 1
            if used_base[b_key] == 0:  # Only write as F if it hasn't been a T
                ret.append(to_process)
            used_base[b_key] = 9
        # Write both (any state)
        elif not base_is_used and not comp_is_used:
            to_process = copy.copy(match)
            # Don't write twice
            if used_base[b_key] != 0:
                to_process.base = None
            if used_comp[c_key] != 0:
                to_process.comp = None

            used_base[b_key] += match.comp_gt_count
            used_comp[c_key] += match.base_gt_count
            # All used up
            if used_base[b_key] >= match.base_gt_count:
                base_cnt -= 1
            if used_comp[c_key] >= match.comp_gt_count:
                comp_cnt -= 1

            # Safety edge case check
            if to_process.base is not None or to_process.comp is not None:
                ret.append(to_process)
    return ret


def pick_single_matches(match_matrix):
    """
    Given a numpy array of MatchResults, find the single best match for calls
    Once all best pairs yielded, the unpaired calls are set to FP/FN and yielded
    """
    hit_order = np.array(np.unravel_index(np.argsort(
        match_matrix, axis=None), match_matrix.shape)).T
    top_hit_idx = len(hit_order) - 1

    mask = np.ones(match_matrix.shape, dtype='bool')  # True = can use
    ret = []

    used_pairs = []
    while top_hit_idx >= 0 and mask.any():
        idx = tuple(hit_order[top_hit_idx])
        if mask[idx]:
            mask[idx[0], :] = False
            mask[:, idx[1]] = False
            used_pairs.append(idx)
            ret.append(match_matrix[idx])
        top_hit_idx -= 1

    used_base, used_comp = zip(*used_pairs) if used_pairs else ([], [])
    # FNs
    for base_col in set(range(match_matrix.shape[0])) - set(used_base):
        comp_col = match_matrix[base_col].argmax()
        to_process = copy.copy(match_matrix[base_col, comp_col])
        to_process.comp = None  # The comp will be written elsewhere
        to_process.multi = to_process.state
        to_process.state = False
        ret.append(to_process)

    # FPs
    for comp_col in set(range(match_matrix.shape[1])) - set(used_comp):
        base_col = match_matrix[:, comp_col].argmax()
        to_process = copy.copy(match_matrix[base_col, comp_col])
        to_process.base = None  # The base will be written elsewhere
        to_process.multi = to_process.state
        to_process.state = False
        ret.append(to_process)
    return ret


PICKERS = {"single": pick_single_matches,
           "ac": pick_ac_matches,
           "multi": pick_multi_matches
           }


def bench_main(cmdargs):
    """
    Main - entry point from command line
    """
    args = parse_args(cmdargs)

    if check_params(args) or check_inputs(args):
        sys.stderr.write("Couldn't run Truvari. Please fix parameters\n")
        sys.exit(100)

    matcher = truvari.Matcher(args)

    m_bench = Bench(matcher, args.base, args.comp, args.output,
                    args.includebed, args.extend, args.debug, True)
    output = m_bench.run()

    logging.info("Stats: %s", json.dumps(output.stats_box, indent=4))
    logging.info("Finished bench")
