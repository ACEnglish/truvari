#!/usr/bin/env python
from __future__ import print_function
import os
import re
import sys
import json
import bisect
import logging
import argparse
import warnings

from collections import defaultdict, Counter

# External dependencies
import vcf
import pyfaidx
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


def make_interval_tree(vcf_file, sizemin=10, sizemax=100000, passonly=False):
    """
    Return a dictonary of {chr:[start, end,..], ...}
    that can be queried with bisect to return a span to fetch variants against
    could stand to make a progress bar since this now takes so long
    """
    n_entries = 0
    cmp_entries = 0
    lookup = defaultdict(IntervalTree)
    try:
        for entry in vcf_file:
            n_entries += 1
            if passonly and (entry.FILTER is not None and len(entry.FILTER)):
                continue
            start, end = get_vcf_boundaries(entry)
            sz = get_vcf_entry_size(entry)
            if sz < sizemin or sz > sizemax:
                continue
            cmp_entries += 1
            lookup[entry.CHROM].addi(start, end, entry.start)
    except ValueError as e:
        logging.error("Unable to parse comparison vcf file. Please check header definitions")
        logging.error("Specific error: \"%s\"", str(e))
        exit(1)
    return lookup, n_entries, cmp_entries


def vcf_to_key(source,  entry):
    """
    Turn a vcf entry into a hashable key using the 'source' (base/comp) to separate the two
    setsource.chr:pos:ref:alt
    helpful for not re-using variants
    BUG: if a caller redundantly calls a variant exactly the same. It will be collapsed
    the
    """
    start, end = get_vcf_boundaries(entry)
    return "%s.%s:%d-%d(%s|%s)" % (source, entry.CHROM, start, end, entry.REF, str(entry.ALT[0]))


def var_sizesim(sizeA, sizeB):
    """
    Calculate the size similarity pct for the two entries
    compares the longer of entryA's two alleles (REF or ALT)
    """
    return min(sizeA, sizeB) / float(max(sizeA, sizeB)), sizeA - sizeB


def gt_comp(entryA, entryB, sampleA, sampleB):
    """
    Compare the genotypes, returns if they're the same
    Simple for now. Methodized so it's easy to expand later
    """
    return entryA.genotype(sampleA)["GT"] == entryB.genotype(sampleB)["GT"]



def create_haplotype(entryA, entryB, ref):
    """
    Turn two entries into their haplotype sequence for comparison
    """
    if len(entryA.ALT[0]) >= len(entryA.REF):
        a1_chrom, a1_start, a1_end, a1_seq = entryA.CHROM, entryA.start, entryA.end, str(entryA.ALT[0]).upper()
        a2_chrom, a2_start, a2_end, a2_seq = entryB.CHROM, entryB.start, entryB.end, str(entryB.ALT[0]).upper()
    else:
        a1_chrom, a1_start, a1_end, a1_seq = entryA.CHROM, entryA.start, entryA.end, entryA.REF  # str(entryA.ALT[0]).upper()
        a2_chrom, a2_start, a2_end, a2_seq = entryB.CHROM, entryB.start, entryB.end, entryB.REF  # str(entryB.ALT[0]).upper()

    start = min(a1_start, a2_start)
    end = max(a1_end, a2_end)

    hap1_seq = ref.get_seq(a1_chrom, start + 1, a1_start).seq + a1_seq + ref.get_seq(a1_chrom, a1_end + 1, end).seq
    hap2_seq = ref.get_seq(a2_chrom, start + 1, a2_start).seq + a2_seq + ref.get_seq(a2_chrom, a2_end + 1, end).seq
    return str(hap1_seq), str(hap2_seq)


def var_pctsim_lev(entryA, entryB, ref):
    """
    Use Levenshtein distance ratio of the larger sequence as a proxy
    to pct sequence similarity
    """
    # Shortcut to save compute - probably unneeded
    if entryA.REF == entryB.REF and entryA.ALT[0] == entryB.ALT[0]:
        return 1.0
    # Handling of breakends should be here
    if entryA.var_subtype == "complex" and entryB.var_subtype == "complex":
        return 1
    elif entryA.var_subtype == "complex" or entryB.var_subtype == "complex":
        return 0
    allele1, allele2 = create_haplotype(entryA, entryB, ref)
    return Levenshtein.ratio(allele1, allele2)


def overlaps(s1, e1, s2, e2):
    """
    Check if two ranges have overlap
    """
    s_cand = max(s1, s2)
    e_cand = min(e1, e2)
    return s_cand < e_cand


def same_variant_type(entryA, entryB):
    """
    Look at INFO/SVTYPE to see if they two calls are of the same type.
    If SVTYPE is unavailable, Infer if this is a insertion or deletion by
    looking at the REF/ALT sequence size differences
    If REF/ALT are not available, try to use the <INS> <DEL> in the ALT column.
    else rely on vcf.model._Record.var_subtype

    return a_type == b_type
    """
    sv_alt_match = re.compile("\<(?P<SVTYPE>.*)\>")

    def pull_type(entry):
        ret_type = None
        if "SVTYPE" in entry.INFO:
            ret_type = entry.INFO["SVTYPE"]
            if type(ret_type) is list:
                logging.warning("SVTYPE is list for entry %s", str(entry))
                ret_type = ret_type[0]
            return ret_type
        if not (str(entry.ALT[0]).count("<" or str(entry.ALT)[0].count(":"))):
            # Doesn't have <INS> or BNDs as the alt seq, then we can assume it's sequence resolved..?
            if len(entry.REF) <= len(entry.ALT[0]):
                ret_type = "INS"
            elif len(entry.REF) >= len(entry.ALT[0]):
                ret_type = "DEL"
            elif len(entry.REF) == len(entry.ALT[0]):
                # Is it really?
                ret_type = "COMPLEX"
            return ret_type
        mat1 = sv_alt_match.match(entry.ALT[0])
        if mat is not None:
            return mat1.groups()["SVTYPE"]
        # rely on pyvcf
        return entry.var_subtype.upper()

    a_type = pull_type(entryA)
    b_type = pull_type(entryB)
    return a_type == b_type


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
        if type(end) == list:
            end = end[0]
    else:
        end = entry.end
    if start > end:
        start, end = end, start
    return start, end


def get_vcf_entry_size(entry):
    """
    Calculate the size of the variant. Use SVLEN INFO tag if available. Otherwise infer
    Return absolute value of the size
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
    return (2 * sim + 1 * size + 1 * ovl) / 3.0


def setup_progressbar(size):
    """
    Return a formatted progress bar of size
    """
    return progressbar.ProgressBar(redirect_stdout=True, max_value=size, widgets=[
        ' [', progressbar.Timer(), ' ', progressbar.Counter(), '/', str(size), '] ',
        progressbar.Bar(),
        ' (', progressbar.ETA(), ') ',
    ])


def annotate_tp(entry, score, pctsim, pctsize, pctovl, szdiff, stdist, endist, oentry):
    """
    Add the matching annotations to a vcf entry
    match_score, match_pctsim, match_pctsize, match_ovlpct, match_szdiff, \
                    match_stdist, match_endist, match_entry
    """
    entry.INFO["PctSeqSimilarity"] = pctsim
    entry.INFO["PctSizeSimilarity"] = pctsize

    entry.INFO["PctRecOverlap"] = pctovl

    entry.INFO["SizeDiff"] = szdiff

    entry.INFO["StartDistance"] = stdist
    entry.INFO["EndDistance"] = endist


def make_giabreport(args, stats_box):
    """
    Create summaries of the TPs/FNs
    """
    def make_entries(vcf_file):
        """
        Turn all vcf entries into a list of dicts for easy parsing
        """
        ret = []
        vcf_reader = vcf.Reader(filename=vcf_file)
        for entry in vcf_reader:
            data = {}
            for key in entry.INFO:
                data[key] = entry.INFO[key]
            for sample in entry.samples:
                for fmt in sample.data._fields:
                    data[sample.sample + "_" + fmt] = sample[fmt]
            data["CHROM"] = entry.CHROM
            data["POS"] = entry.POS
            data["ID"] = entry.ID
            if entry.FILTER is not None:
                data["FILTER"] = ";".join(entry.FILTER)
            data["QUAL"] = entry.QUAL
            ret.append(data)
        return ret

    def count_by(key, docs, out):
        """
        get a count for all the keys
        key is a dictionary of count and their order
        """
        main_key = key.keys()[0]
        cnt = Counter()
        for x in docs:
            cnt[x[main_key]] += 1
        for k in key[main_key]:
            out.write("%s\t%s\n" % (k, cnt[k]))

    def twoxtable(key1, key2, docs, out):
        """
        Parse any set of docs and create a 2x2 table
        """
        main_key1 = key1.keys()[0]
        main_key2 = key2.keys()[0]
        cnt = defaultdict(Counter)
        for x in docs:
            cnt[x[main_key1]][x[main_key2]] += 1

        out.write(".\t" + "\t".join([str(i) for i in key1[main_key1]]) + '\n')
        for y in key2[main_key2]:
            o = [str(y)]
            for x in key1[main_key1]:
                o.append(str(cnt[x][y]))
            out.write("\t".join(o) + '\n')

    def bool_counter(keys, docs, out):
        """
        For bool valued INFO keys, write summaries of their counts
        """
        cnt = defaultdict(Counter)
        col_names = set()
        for x in docs:
            for y in keys:
                cnt[y][x[y]] += 1
                col_names.add(x[y])
        col_names = sorted(list(col_names))
        out.write(".\t%s\n" % "\t".join(col_names))
        row_keys = sorted(list(cnt.keys()))
        for row in row_keys:
            out.write(row)
            for col_val in col_names:
                out.write("\t%d" % cnt[row][col_val])
            out.write("\n")

    def collapse_techs(docs):
        """
        Make a new annotation about the presence of techs inplace
        called techs
        Illcalls
            PBcalls
        CGcalls
        TenXcalls
        """
        calls = ["Illcalls", "PBcalls", "CGcalls", "TenXcalls"]
        for d in docs:
            new_anno = []
            for i in calls:
                if d[i] > 0:
                    new_anno.append(i.rstrip("calls"))
            d["techs"] = "+".join(new_anno)

    logging.info("Creating GIAB report")

    sum_out = open(os.path.join(args.output, "giab_report.txt"), 'w')

    tp_base = make_entries(os.path.join(args.output, "tp-base.vcf"))

    collapse_techs(tp_base)
    fn = make_entries(os.path.join(args.output, "fn.vcf"))
    collapse_techs(fn)

    size_keys = {"sizecat": ["50to99", "100to299", "300to999", "gt1000"]}
    svtype_keys = {"SVTYPE": ["DEL", "INS", "COMPLEX"]}
    tech_keys = {"techs": ["I+PB+CG+TenX", "I+PB+CG", "I+PB+TenX", "PB+CG+TenX",
                           "I+PB", "I+CG", "I+TenX", "PB+CG", "PB+TenX", "CG+TenX",
                           "I", "PB", "CG", "TenX"]}
    rep_keys = {"REPTYPE": ["SIMPLEDEL", "SIMPLEINS", "DUP", "SUBSDEL", "SUBSINS", "CONTRAC"]}
    tr_keys = ["TRall", "TRgt100", "TRgt10k", "segdup"]

    gt_keys_proband = {"HG002_GT": ["0/1", "./1", "1/1"]}
    gt_keys_father = {"HG003_GT": ["./.", "0/0", "0/1", "./1", "1/1"]}
    gt_keys_mother = {"HG004_GT": ["./.", "0/0", "0/1", "./1", "1/1"]}
    # OverallNumbers
    sum_out.write("TP\t%s\n" % (len(tp_base)))
    sum_out.write("FN\t%s\n\n" % (len(fn)))
    sum_out.write("TP_size\n")
    count_by(size_keys, tp_base, sum_out)
    sum_out.write("FN_size\n")
    count_by(size_keys, fn, sum_out)
    sum_out.write("\nTP_type\n")
    count_by(svtype_keys, tp_base, sum_out)
    sum_out.write("FN_type\n")
    count_by(svtype_keys, fn, sum_out)
    sum_out.write("\nTP_Type+Size\n")
    twoxtable(svtype_keys, size_keys, tp_base, sum_out)
    sum_out.write("FN_Type+Size\n")
    twoxtable(svtype_keys, size_keys, fn, sum_out)
    sum_out.write("\nTP_REPTYPE\n")
    count_by(rep_keys, tp_base, sum_out)
    sum_out.write("FN_REPTYPE\n")
    count_by(rep_keys, fn, sum_out)
    sum_out.write("\nTP_size+REPTYPE\n")
    twoxtable(size_keys, rep_keys, tp_base, sum_out)
    sum_out.write("FN_size+REPTYPE\n")
    twoxtable(size_keys, rep_keys, fn, sum_out)
    sum_out.write("\nTP_Tech\n")
    count_by(tech_keys, tp_base, sum_out)
    sum_out.write("FN_Tech\n")
    count_by(tech_keys, fn, sum_out)
    sum_out.write("\nTP_Size+Tech\n")
    twoxtable(size_keys, tech_keys, tp_base, sum_out)
    sum_out.write("FN_Size+Tech\n")
    twoxtable(size_keys, tech_keys, fn, sum_out)
    sum_out.write("\nTP_Type+Tech\n")
    twoxtable(svtype_keys, tech_keys, tp_base, sum_out)
    sum_out.write("FN_Type+Tech\n")
    twoxtable(svtype_keys, tech_keys, fn, sum_out)

    sum_out.write("\nPerformance\n")
    # Add output of the stats box
    for key in sorted(stats_box.keys()):
        sum_out.write("%s\t%s\n" % (key, str(stats_box[key])))

    sum_out.write("\nArgs\n")
    # Add output of the parameters
    argd = vars(args)
    for key in sorted(argd.keys()):
        sum_out.write("%s\t%s\n" % (key, str(argd[key])))
    #
    sum_out.write("\nTP_HG002GT\n")
    count_by(gt_keys_proband, tp_base, sum_out)
    sum_out.write("FN_HG002GT\n")
    count_by(gt_keys_proband, fn, sum_out)

    sum_out.write("\nTP_HG003.HG004GT\n")
    twoxtable(gt_keys_father, gt_keys_mother, tp_base, sum_out)
    sum_out.write("FN_HG003.HG004GT\n")
    twoxtable(gt_keys_father, gt_keys_mother, fn, sum_out)

    sum_out.write("\nTP TandemRepeat Anno\n")
    bool_counter(tr_keys, tp_base, sum_out)
    sum_out.write("FN TandemRepeat Anno\n")
    bool_counter(tr_keys, fn, sum_out)
    sum_out.close()


class GenomeTree():

    """
    Helper class to specify included regions of the genome when iterating events.
    """

    def __init__(self, vcfA, vcfB, include=None):
        contigA_set = set(vcfA.contigs.keys())
        contigB_set = set(vcfB.contigs.keys())
        all_regions = defaultdict(IntervalTree)
        if include is not None:
            counter = 0
            with open(include, 'r') as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    data = line.strip().split('\t')
                    chrom = data[0]
                    start = int(data[1])
                    end = int(data[2])
                    all_regions[chrom].addi(start, end)
                    counter += 1
            logging.info("Including %d bed regions", counter)
        else:
            excluding = contigB_set - contigA_set
            if len(excluding):
                logging.warning(
                    "Excluding %d contigs present in comparison calls header but not base calls.", len(excluding))

            for contig in contigA_set:
                name = vcfA.contigs[contig].id
                length = vcfA.contigs[contig].length
                all_regions[name].addi(0, length)

        self.tree = all_regions

    def iterate(self, vcf_file):
        """
        Iterates a vcf and yields only the entries that overlap an 'include' region
        """
        for entry in vcf_file:
            if self.include(entry):
                yield entry

    def include(self, entry):
        """
        Returns if this entry's start and end are within a region that is to be included
        Here overlap means lies completely within the boundary of an include region
        """
        astart, aend = get_vcf_boundaries(entry)
        overlaps = self.tree[entry.CHROM].overlaps(astart) and self.tree[entry.CHROM].overlaps(aend)
        if astart == aend:
            return overlaps
        return overlaps and len(self.tree[entry.CHROM].search(astart, aend)) == 1


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
    my_vcf.infos['PctRecOverlap'] = vcf.parser._Info(
        id="PctRecOverlap", num=1, type='Float',
        desc="Percent reciprocal overlap percent of the two calls' coordinates",
        source=None, version=None)
    my_vcf.infos['StartDistance'] = vcf.parser._Info(
        id="StartDistance", num=1, type='Integer',
        desc="Distance of this call's start from comparison call's start",
        source=None, version=None)
    my_vcf.infos['EndDistance'] = vcf.parser._Info(
        id="EndDistance", num=1, type='Integer',
        desc="Distance of this call's start from comparison call's start",
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
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
        return x

    parser = argparse.ArgumentParser(prog="truvari", description=USAGE,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--base", type=str, required=True,
                        help="Baseline truth-set calls")
    parser.add_argument("-c", "--comp", type=str, required=True,
                        help="Comparison set of calls")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output directory")
    parser.add_argument("-f", "--reference", type=str, required=True, default=None,
                        help="Indexed fasta used to call variants")
    parser.add_argument("--giabreport", action="store_true",
                        help="Parse output TPs/FNs for GIAB annotations and create a report")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Verbose logging")

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
    filteg.add_argument("--no-ref", action="store_true", default=False,
                        help="Don't include 0/0 or ./. GT calls (%(default)s)")
    filteg.add_argument("--includebed", type=str, default=None,
                        help="Bed file of regions in the genome to include only calls overlapping")
    filteg.add_argument("--multimatch", action="store_true", default=False,
                        help="Allow base calls to match multiple comparison calls (%(default)s)")

    args = parser.parse_args(args)
    return args


def run(cmdargs):
    args = parse_args(cmdargs)
    if os.path.isdir(args.output):
        print("Error! Output directory '%s' already exists" % args.output)
        exit(1)

    reference = pyfaidx.Fasta(args.reference)
    os.mkdir(args.output)

    setup_logging(args.debug, LogFileStderr(os.path.join(args.output, "log.txt")))
    logging.info("Params:\n%s", json.dumps(vars(args), indent=4))

    ref_gts = ["0/0", "0|0", "./.", ".|."]

    vcf_base = vcf.Reader(filename=args.base)
    if args.bSample is not None:
        sampleA = args.bSample
        if sampleA not in vcf_base.samples:
            logging.error("Sample %s not found in vcf (%s)", sampleA, vcf_base.samples)
            exit(1)
    else:
        sampleA = vcf_base.samples[0]

    # Check early so we don't waste users' time
    if not os.path.exists(args.comp):
        logging.error("File %s does not exist" % (args.comp))
        exit(1)
    if not os.path.exists(args.base):
        logging.error("File %s does not exist" % (args.base))
        exit(1)
    check_fail = False
    if not args.comp.endswith(".gz"):
        check_fail = True
        logging.error("Comparison vcf %s does not end with .gz. Must be bgzip'd", args.comp)
    if not os.path.exists(args.comp + '.tbi'):
        check_fail = True
        logging.error("Comparison vcf index %s.tbi does not exist. Must be indexed", args.comp)
    if check_fail:
        exit(1)

    vcf_comp = vcf.Reader(filename=args.comp)
    edit_header(vcf_comp)
    if args.cSample is not None:
        sampleB = args.cSample
        if sampleB not in vcf_comp.samples:
            logging.error("Sample %s not found in vcf (%s)", sampleB, vcf_comp.samples)
            exit(1)
    else:
        sampleB = vcf_comp.samples[0]

    regions = GenomeTree(vcf_base, vcf_comp, args.includebed)

    logging.info("Creating call interval tree for overlap search")
    span_lookup, num_entries_b, cmp_entries = make_interval_tree(vcf_comp, args.sizefilt, args.sizemax, args.passonly)

    logging.info("%d call variants", num_entries_b)
    logging.info("%d call variants within size range (%d, %d)", cmp_entries, args.sizefilt, args.sizemax)

    num_entries = 0
    for entry in regions.iterate(vcf_base):
        num_entries += 1
    logging.info("%s base variants", num_entries)
    # Reset
    vcf_base = vcf.Reader(filename=args.base)
    edit_header(vcf_base)

    bar = setup_progressbar(num_entries)

    # Setup outputs
    tpb_out = vcf.Writer(open(os.path.join(args.output, "tp-base.vcf"), 'w'), vcf_base)
    b_filt = vcf.Writer(open(os.path.join(args.output, "base-filter.vcf"), 'w'), vcf_base)
    tpc_out = vcf.Writer(open(os.path.join(args.output, "tp-call.vcf"), 'w'), vcf_comp)
    c_filt = vcf.Writer(open(os.path.join(args.output, "call-filter.vcf"), 'w'), vcf_comp)
    fn_out = vcf.Writer(open(os.path.join(args.output, "fn.vcf"), 'w'), vcf_base)
    fp_out = vcf.Writer(open(os.path.join(args.output, "fp.vcf"), 'w'), vcf_comp)

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

    # Calls that have been matched up
    matched_calls = defaultdict(bool)
    # for variant A - do var_match in B
    logging.info("Matching base to calls")
    for pbarcnt, base_entry in enumerate(regions.iterate(vcf_base)):
        bar.update(pbarcnt)
        sizeA = get_vcf_entry_size(base_entry)

        if sizeA < args.sizemin or sizeA > args.sizemax:
            stats_box["base size filtered"] += 1
            b_filt.write_record(base_entry)
            continue
        if args.no_ref and not base_entry.genotype(sampleA).is_variant:
            stats_box["base gt filtered"] += 1
            b_filt.write_record(base_entry)
            continue

        if args.passonly and (base_entry.FILTER is not None and len(base_entry.FILTER)):
            continue

        stats_box["base cnt"] += 1

        fetch_start, fetch_end = fetch_coords(span_lookup, base_entry, args.refdist)
        if fetch_start is None and fetch_end is None:
            # No overlaps, don't even bother checking
            base_entry.INFO["NumNeighbors"] = 0
            base_entry.INFO["NumThresholdNeighbors"] = 0
            stats_box["FN"] += 1
            fn_out.write_record(base_entry)
            continue

        # IntervalTree can give boundaries past REFDIST in the case of Inversions where start>end
        # We still need to fetch on the expanded boundaries so we can test them, but
        #   we need to filter calls that otherwise shouldn't be considered
        #  see the bstart/bend below
        astart, aend = get_vcf_boundaries(base_entry)

        thresh_neighbors = []
        num_neighbors = 0

        # methodize so we can parallelize?
        # - i'm woried the vcfRecord isn't pickleable so we'd need to collect them and then return
        # more importantly, we'd also have race conditions for when it's been used...
        # even if already_considered was atomic, you may get diff answers
        # +- 1 just to be safe because why not...
        for comp_entry in vcf_comp.fetch(base_entry.CHROM, max(0, fetch_start - 1), fetch_end + 1):

            # There is a race condition here that could potentially mismatch things
            # If base1 passes matching call1 and then base2 passes matching call1
            # better, it can't use it and we mismatch -- UPDATE: by default we don't enforce one-match
            logging.debug("Comparing %s %s", str(base_entry), str(comp_entry))
            b_key = vcf_to_key('c', comp_entry)
            if not args.multimatch and matched_calls[b_key]:
                logging.debug("No match because comparison call already matched")
                continue

            sizeB = get_vcf_entry_size(comp_entry)
            if sizeB < args.sizefilt:
                continue

            # Double ensure OVERLAP - there's a weird edge case where fetch with
            # the interval tree can return non-overlaps
            bstart, bend = get_vcf_boundaries(comp_entry)
            if not overlaps(astart - args.refdist, aend + args.refdist, bstart, bend):
                continue
            if not regions.include(comp_entry):
                continue
            # Someone in the Base call's neighborhood, we'll see if it passes comparisons
            num_neighbors += 1

            if args.no_ref and not comp_entry.genotype(sampleB).is_variant:
                logging.debug("%s is hom-ref", str(comp_entry))
                continue

            if args.gtcomp and not gt_comp(base_entry, comp_entry, sampleA, sampleB):
                logging.debug("%s and %s are not the same genotype", str(base_entry), str(comp_entry))
                continue

            if not args.typeignore and not same_variant_type(base_entry, comp_entry):
                logging.debug("%s and %s are not the same SVTYPE", str(base_entry), str(comp_entry))
                continue

            size_similarity, size_diff = var_sizesim(sizeA, sizeB)
            if size_similarity < args.pctsize:
                logging.debug("%s and %s size similarity is too low (%f)", str(
                    base_entry), str(comp_entry), size_similarity)
                continue

            ovl_pct = get_rec_ovl(astart, aend, bstart, bend)
            if ovl_pct < args.pctovl:
                logging.debug("%s and %s overlap percent is too low (%f)", str(base_entry), str(comp_entry), ovl_pct)
                continue

            if args.pctsim > 0:
                seq_similarity = var_pctsim_lev(base_entry, comp_entry, reference)
                if seq_similarity < args.pctsim:
                    logging.debug("%s and %s sequence similarity is too low (%f)", str(
                        base_entry), str(comp_entry), seq_similarity)
                    continue
            else:
                seq_similarity = 0

            start_distance = astart - bstart
            end_distance = aend - bend

            score = get_weighted_score(seq_similarity, size_similarity, ovl_pct)
            # If you put these numbers in an object, it'd be easier to pass round
            # You'd just need to make it sortable
            thresh_neighbors.append(
                (score, seq_similarity, size_similarity, ovl_pct, size_diff, start_distance, end_distance, comp_entry))

        # reporting the best match
        base_entry.INFO["NumNeighbors"] = num_neighbors
        base_entry.INFO["NumThresholdNeighbors"] = len(thresh_neighbors)
        if len(thresh_neighbors) > 0:
            thresh_neighbors.sort(reverse=True)
            logging.debug("Picking from candidate matches:\n%s", "\n".join([str(x) for x in thresh_neighbors]))

            match_score, match_pctsim, match_pctsize, match_ovlpct, match_szdiff, \
                match_stdist, match_endist, match_entry = thresh_neighbors[0]
            logging.debug("Best match is %s", str(match_entry))

            base_entry.INFO["TruScore"] = match_score
            match_entry.INFO["TruScore"] = match_score
            match_entry.INFO["NumNeighbors"] = num_neighbors
            match_entry.INFO["NumThresholdNeighbors"] = len(thresh_neighbors)

            # Don't double count calls was found before
            b_key = vcf_to_key('b', base_entry)
            c_key = vcf_to_key('c', match_entry)
            if not matched_calls[b_key]:
                stats_box["TP-base"] += 1
            if not matched_calls[c_key]:
                stats_box["TP-call"] += 1
                write_tp_call = True
            else:
                write_tp_call = False
            # Mark the call for multimatch checking
            matched_calls[b_key] = True
            matched_calls[c_key] = True

            annotate_tp(base_entry, *thresh_neighbors[0])
            annotate_tp(match_entry, *thresh_neighbors[0])

            tpb_out.write_record(base_entry)
            if write_tp_call:
                tpc_out.write_record(match_entry)
            logging.debug("Matching %s and %s", str(base_entry), str(match_entry))
        else:
            stats_box["FN"] += 1
            fn_out.write_record(base_entry)

    bar.finish()

    # Get a results peek
    do_stats_math = True
    if stats_box["TP-base"] == 0 and stats_box["FN"] == 0:
        logging.warning("No TP or FN calls in base!")
        do_stats_math = False
    else:
        logging.info("Results peek: %d TP-base %d FN %.2f%% Recall", stats_box["TP-base"], stats_box["FN"],
                     100 * (float(stats_box["TP-base"]) / (stats_box["TP-base"] + stats_box["FN"])))

    logging.info("Parsing FPs from calls")
    bar = setup_progressbar(num_entries_b)
    # Reset
    vcf_comp = vcf.Reader(filename=args.comp)
    edit_header(vcf_comp)
    for cnt, entry in enumerate(regions.iterate(vcf_comp)):
        # Need to count these, I think
        if args.passonly and (entry.FILTER is not None and len(entry.FILTER)):
            continue
        bar.update(cnt + 1)
        if matched_calls[vcf_to_key('c', entry)]:
            continue
        size = get_vcf_entry_size(entry)
        if size < args.sizemin or size > args.sizemax:
            c_filt.write_record(entry)
            stats_box["call size filtered"] += 1
        elif args.no_ref and not entry.genotype(sampleB).is_variant:
            stats_box["call gt filtered"] += 1
        elif regions.include(entry):
            fp_out.write_record(entry)
            stats_box["FP"] += 1
    bar.finish()
    # call count is just of those used were used
    stats_box["call cnt"] = stats_box["TP-base"] + stats_box["FP"]

    # Close to flush vcfs
    tpb_out.close()
    b_filt.close()
    tpc_out.close()
    c_filt.close()
    fn_out.close()
    fp_out.close()

    # Final calculations
    if do_stats_math:
        # precision
        stats_box["precision"] = float(stats_box["TP-call"]) / (stats_box["TP-call"] + stats_box["FP"])
        # recall
        stats_box["recall"] = float(stats_box["TP-base"]) / (stats_box["TP-base"] + stats_box["FN"])
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
        sys.stdout.write(json.dumps(stats_box, indent=4) + '\n')
        logging.info("Stats: %s", json.dumps(stats_box, indent=4))

    if args.giabreport:
        make_giabreport(args, stats_box)

    logging.info("Finished")

if __name__ == '__main__':
    run(sys.argv[1:])
