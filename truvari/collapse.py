"""
Structural variant collapser

Will collapse all variants within sizemin/max that match over thresholds
"""
import os
import sys
import gzip
import json
import logging
import argparse
import itertools
import statistics
from dataclasses import dataclass, field
from functools import cmp_to_key, partial

import numpy as np

import truvari
import truvari.bench as trubench


@dataclass
class CollapsedCalls():
    """
    Holds all relevant information for a set of collapsed calls
    """
    entry: None
    match_id: str
    matches: list = field(default_factory=list)
    gt_consolidate_count: int = 0
    genotype_mask: str = None  # not actually a str

    def calc_median_sizepos(self):
        """
        Given a set of variants, calculate the median start, end, and size
        """
        st, en = self.entry.boundaries()
        sz = self.entry.var_size()
        starts = [st]
        ends = [en]
        sizes = [sz]
        for i in self.matches:
            st, en = i.comp.boundaries()
            sz = i.comp.var_size()
            starts.append(st)
            ends.append(en)
            sizes.append(sz)
        return int(statistics.median(starts)), int(statistics.median(ends)), int(statistics.median(sizes))

    def annotate_entry(self, header, med_info):
        """
        Edit an entry that's going to be collapsed into another entry
        """
        self.entry.translate(header)
        self.entry.info["NumCollapsed"] = len(self.matches)
        self.entry.info["NumConsolidated"] = self.gt_consolidate_count
        self.entry.info["CollapseId"] = self.match_id
        if med_info:
            self.entry.info["CollapseStart"], self.entry.info["CollapseEnd"], self.entry.info["CollapseSize"] = self.calc_median_sizepos()

    @staticmethod
    def make_genotype_mask(entry, gtmode):
        """
        Populate the genotype mask
        """
        if gtmode == 'off':
            return None
        to_mask = (lambda x: 1 in x) if gtmode == 'all' else (
            lambda x: x.count(1) == 1)
        return np.array([to_mask(_.allele_indices) for _ in entry.samples.values()], dtype=bool)

    def gt_conflict(self, other, which_gt):
        """
        Return true if entry's genotypes conflict with any of the current collapse
        which_gt all prevents variants present in the same sample from being collapsed
        which_gt het only prevents two het variants from being collapsed.
        """
        if which_gt == 'off':
            return False

        o_mask = self.make_genotype_mask(other, which_gt)
        if (self.genotype_mask & o_mask).any():
            return True
        self.genotype_mask |= o_mask
        return False


def find_new_matches(base, remaining_calls, dest_collapse, params):
    """
    Pull variants from remaining calls that match to the base entry into the destination collapse
    Updates everything in place
    """
    # Sort based on size difference to current call
    remaining_calls.sort(key=partial(relative_size_sorter, base))
    i = 0
    while i < len(remaining_calls):
        candidate = remaining_calls[i]
        mat = base.match(candidate)
        mat.matid = dest_collapse.match_id
        if params.hap and not hap_resolve(base,  candidate):
            mat.state = False
        if mat.state and dest_collapse.gt_conflict(candidate, params.gt):
            mat.state = False
        if mat.state:
            dest_collapse.matches.append(mat)
            remaining_calls.pop(i)
        else:
            # move to next one
            i += 1
            # short circuit
            if mat.sizesim is not None and mat.sizesim < params.pctsize:
                return

def collapse_chunk(chunk, params):
    """
    Returns a list of lists with [keep entry, collap matches, match_id, gt_consolidate_count]
    """
    chunk_dict, chunk_id = chunk
    remaining_calls = sorted(chunk_dict['base'], key=params.sorter)

    call_id = -1
    ret = []  # list of Collapses
    while remaining_calls:
        call_id += 1
        m_collap = CollapsedCalls(remaining_calls.pop(0),
                                  f'{chunk_id}.{call_id}')
        # quicker genotype comparison - needs to be refactored
        m_collap.genotype_mask = m_collap.make_genotype_mask(m_collap.entry,
                                                             params.gt)

        find_new_matches(m_collap.entry, remaining_calls, m_collap, params)
        # Chain will only operate once to prevent 'sliding'
        # e.g. collapsing integers within 5 of 1, 4, 6, 8
        # without a single operation 1
        # with a single operation 1, 8
        # not chaining at all 1, 6
        if params.chain:
            i = 0
            while remaining_calls and i < len(m_collap.matches):
                find_new_matches(m_collap.matches[i].comp, remaining_calls, m_collap, params)
                i += 1

        # If hap, only allow the best match
        # put the rest of the remaining calls back
        if params.hap and m_collap.matches:
            mats = sorted(m_collap.matches, reverse=True)
            m_collap.matches = [mats.pop(0)]
            remaining_calls.extend(mats)
        ret.append(m_collap)
        # Sort back to where they need to be to choose the next to evaluate
        remaining_calls.sort(key=params.sorter)

    if params.no_consolidate:
        for val in ret:
            if params.gt != 'off':
                edited_entry, collapse_cnt = gt_aware_consolidate(
                    val.entry, val.matches)
            else:
                edited_entry, collapse_cnt = collapse_into_entry(
                    val.entry, val.matches, params.hap)
            val.entry = edited_entry
            val.gt_consolidate_count = collapse_cnt

    for i in chunk_dict['__filtered']:
        ret.append(CollapsedCalls(i, None))
    ret.sort(key=cmp_to_key(lambda x, y: x.entry.pos - y.entry.pos))
    return ret


def relative_size_sorter(base, comp):
    """
    Sort calls based on the absolute size difference of base and comp
    """
    return abs(base.var_size() - comp.var_size())


def collapse_into_entry(entry, others, hap_mode=False):
    """
    Consolidate information for genotypes where sample is unset
    okay - this works, but its a mess
    """
    # short circuit
    if not others:
        return entry, 0

    # We'll populate with the most similar, first
    others.sort(reverse=True)
    # I have a special case of --hap. I need to allow hets
    replace_gts = [truvari.GT.REF, truvari.GT.NON, truvari.GT.UNK]
    if hap_mode:
        replace_gts.insert(1, truvari.GT.HET)

    # Each sample of this entry needs to be checked/set
    n_consolidate = 0
    for sample in entry.samples:
        m_gt = truvari.get_gt(entry.samples[sample]["GT"])
        if m_gt not in replace_gts:
            continue  # already set
        n_idx = None
        o_gt = None
        for pos, o_entry in enumerate(others):
            o_entry = o_entry.comp
            o_gt = truvari.get_gt(o_entry.samples[sample]["GT"])
            if o_gt not in replace_gts:
                n_idx = pos
                break  # this is the first other that's set
        # consolidate
        if hap_mode and m_gt == truvari.GT.HET and o_gt == truvari.GT.HET:
            entry.samples[sample]["GT"] = (1, 1)
            n_consolidate += 1
        elif n_idx is not None:
            n_consolidate += 1
            o_entry = others[n_idx].comp
            for key in set(entry.samples[sample].keys() + o_entry.samples[sample].keys()):
                try:
                    entry.samples[sample][key] = o_entry.samples[sample][key]
                except TypeError:  # pragma: no cover
                    # Happens for things like PL when one is null but its expecting a tuple
                    logging.debug("Unable to set FORMAT %s for sample %s",
                                  key, sample)
                    logging.debug("Kept entry: %s:%d %s",
                                  entry.chrom, entry.pos, entry.id)
                    logging.debug("Colap entry: %s:%d %s",
                                  o_entry.chrom, o_entry.pos, o_entry.id)
                except KeyError:  # pragma: no cover
                    logging.debug("Unshared format %s in sample %s ignored for pair %s:%d %s %s:%d %s",
                                  key, sample, entry.chrom, entry.pos, entry.id, o_entry.chrom,
                                  o_entry.pos, o_entry.id)
            # pass along phase
            entry.samples[sample].phased = o_entry.samples[sample].phased
    return entry, n_consolidate


def get_ac(gt):
    """
    Helper method to get allele count. assumes only 1s as ALT
    """
    return sum(1 for _ in gt if _ == 1)


def get_none(entry, key):
    """
    Make a none_tuple for a format
    """
    cnt = entry.header.formats[key].number
    if cnt in [1, 'A', '.']:
        return None
    if cnt == 'R':
        return (None, None)
    return (None, None, None)


def fmt_none(value):
    """
    Checks all values in a format field for nones
    """
    if isinstance(value, (tuple, list)):
        for i in value:
            if i is not None:
                return value
        return None
    return value


def gt_aware_consolidate(entry, others):
    """
    All formats are consolidated (first one taken)
    And two hets consolidated become hom
    Phase is lost
    """
    if not others:
        return entry, 0
    all_fmts = set(entry.samples[0].keys())
    for i in others:
        all_fmts.update(i.comp.samples[0].keys())
    all_fmts.remove('GT')
    all_fmts = sorted(list(all_fmts))
    n_consolidated = 0
    for sample in entry.samples:
        new_fmts = {"GT": 0}
        # Need to consolidate non-none GT if it isn't present here
        if None in entry.samples[sample]['GT']:
            o_gt = None
            i_phased = False
            i = 0
            while o_gt is None and i < len(others):
                m_fmt = others[i].comp.samples[sample]
                if None not in m_fmt['GT']:
                    i_phased = m_fmt.phased
                    o_gt = m_fmt['GT']
                i += 1
            if o_gt:
                entry.samples[sample]['GT'] = o_gt
                entry.samples[sample].phased = i_phased
        new_fmts['GT'] += get_ac(entry.samples[sample]['GT'])
        # consolidate format fields
        for k in all_fmts:
            new_fmts[k] = None if k not in entry.samples[sample].keys(
            ) else fmt_none(entry.samples[sample][k])
        for o in others:
            o = o.comp
            n_c = get_ac(o.samples[sample]['GT'])
            n_consolidated += 1 if n_c != 0 else 0
            new_fmts['GT'] += n_c
            for k in all_fmts:
                rpl_val = None if k not in o.samples[sample] else fmt_none(
                    o.samples[sample][k])
                # Replace blanks in this sample's new formats
                if rpl_val and new_fmts[k] is None:
                    new_fmts[k] = rpl_val
        if new_fmts['GT'] == 0:
            new_fmts['GT'] = (0, 0)
        elif new_fmts['GT'] == 1:
            new_fmts['GT'] = (0, 1)
        else:
            new_fmts['GT'] = (1, 1)
        for key in all_fmts:
            val = new_fmts[key]
            if val is None:
                val = get_none(entry, key)
            entry.samples[sample][key] = val
    return entry, n_consolidated


def hap_resolve(entryA, entryB):
    """
    Returns true if the calls' genotypes suggest it can be collapsed
    i.e. if either call is HOM, they can't be collapsed.
    If calls are on the same haplotype (1/0 & 1/0), they cannot be collapsed
    """
    gtA = entryA.samples[0]["GT"]
    gtB = entryB.samples[0]["GT"]
    if gtA == (1, 1) or gtB == (1, 1):
        return False
    if gtA == gtB:
        return False
    return True


def sort_length(b1, b2):
    """
    Order entries from longest to shortest SVLEN, ties are by alphanumeric of REF
    """
    s1 = b1.var_size()
    s2 = b2.var_size()
    if s1 < s2:
        return 1
    if s1 > s2:
        return -1
    if b1.ref < b2.ref:
        return 1
    if b1.ref > b2.ref:
        return -1
    return 0


def sort_first(b1, b2):
    """
    Order entries from left-most to right-most POS
    """
    if b1.pos > b2.pos:
        return 1
    if b1.pos < b2.pos:
        return -1
    return sort_length(b1, b2)


def sort_maxqual(b1, b2):
    """
    Order entries from highest to lowest qual
    """
    if b1.qual < b2.qual:
        return 1
    if b1.qual > b2.qual:
        return -1
    return sort_first(b1, b2)


def sort_common(b1, b2):
    """
    Order entries from highest to highest AC
    """
    mac1 = b1.allele_freq_annos()["AC"]
    mac2 = b2.allele_freq_annos()["AC"]
    if mac1 > mac2:
        return 1
    if mac1 < mac2:
        return -1
    return sort_first(b1, b2)


SORTS = {'first': cmp_to_key(sort_first),
         'maxqual': cmp_to_key(sort_maxqual),
         'common': cmp_to_key(sort_common)}


#######
# VCF #
#######
def edit_header(my_vcf, median_info=False):
    """
    Add INFO for new fields to vcf
    """
    header = my_vcf.header.copy()
    header.add_line(('##INFO=<ID=NumCollapsed,Number=1,Type=Integer,'
                     'Description="Number of calls collapsed into this call by truvari">'))
    header.add_line(('##INFO=<ID=CollapseId,Number=1,Type=String,'
                     'Description="Truvari uid to help tie output.vcf and output.collapsed.vcf entries together">'))
    header.add_line(('##INFO=<ID=NumConsolidated,Number=1,Type=Integer,'
                     'Description="Number of samples consolidated into this call by truvari">'))
    if median_info:
        header.add_line(('##INFO=<ID=CollapseStart,Number=1,Type=Integer,'
                        'Description="Median start position of collapsed variants">'))
        header.add_line(('##INFO=<ID=CollapseEnd,Number=1,Type=Integer,'
                        'Description="Median end position of collapsed variants">'))
        header.add_line(('##INFO=<ID=CollapseSize,Number=1,Type=Integer,'
                        'Description="Median size of collapsed variants">'))
    return header


##################
# Args & Outputs #
##################
def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="collapse", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Input variants")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output vcf (stdout)")
    parser.add_argument("-c", "--removed-output", type=str, default="removed.vcf",
                        help="Variants that collapsed into kept variants (%(default)s)")
    parser.add_argument("-f", "--reference", type=str, default=None,
                        help="Indexed fasta used to call variants. Only needed with symbolic variants.")
    parser.add_argument("-w", "--write-resolved", action="store_true",
                        help="Replace resolved symbolic variants with their sequence")
    parser.add_argument("-k", "--keep", choices=["first", "maxqual", "common"], default="first",
                        help="When collapsing calls, which one to keep (%(default)s)")
    parser.add_argument("--bed", type=str, default=None,
                        help="Bed file of regions to analyze")
    parser.add_argument("--gt", type=str, choices=['off', 'all', 'het'], default='off',
                        help="Disallow intra-sample events to collapse for genotypes (%(default)s)")
    parser.add_argument("--intra", action="store_true",
                        help="Intrasample merge to first sample in output (%(default)s)")
    parser.add_argument("--median-info", action="store_true",
                        help="Store median start/end/size of collapsed entries in kept's INFO")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Verbose logging")

    thresg = parser.add_argument_group("Comparison Threshold Arguments")
    thresg.add_argument("-r", "--refdist", type=truvari.restricted_int, default=500,
                        help="Max reference location distance (%(default)s)")
    thresg.add_argument("-p", "--pctseq", type=truvari.restricted_float, default=0.95,
                        help="Min percent sequence similarity. Set to 0 to ignore. (%(default)s)")
    thresg.add_argument("-P", "--pctsize", type=truvari.restricted_float, default=0.95,
                        help="Min pct allele size similarity (minvarsize/maxvarsize) (%(default)s)")
    thresg.add_argument("-O", "--pctovl", type=truvari.restricted_float, default=0.0,
                        help="Min pct reciprocal overlap (%(default)s) for DEL events")
    thresg.add_argument("-t", "--typeignore", action="store_true", default=False,
                        help="Variant types don't need to match to compare (%(default)s)")
    thresg.add_argument("-n", "--no-roll", action="store_false",
                        help="Turn off rolling sequence similarity")
    thresg.add_argument("-m", "--max-resolve", type=truvari.restricted_int, default=25000,
                        help="Maximum size of variant to attempt to sequence resolve ($(default)s)")
    thresg.add_argument("-D", "--decompose", action="store_true",
                        help="Allow decomposition for SV to BND comparison (%(default)s)")
    thresg.add_argument("-d", "--dup-to-ins", action="store_true",
                        help="Assume DUP svtypes are INS (%(default)s)")


    parser.add_argument("--hap", action="store_true", default=False,
                        help="Collapsing a single individual's haplotype resolved calls (%(default)s)")
    parser.add_argument("--chain", action="store_true", default=False,
                        help="Chain comparisons to extend possible collapsing (%(default)s)")
    parser.add_argument("--no-consolidate", action="store_false", default=True,
                        help="Skip consolidation of sample genotype fields (%(default)s)")
    filteg = parser.add_argument_group("Filtering Arguments")
    filteg.add_argument("-s", "--sizemin", type=truvari.restricted_int, default=50,
                        help="Minimum variant size to consider for comparison (%(default)s)")
    filteg.add_argument("-S", "--sizemax", type=truvari.restricted_int, default=50000,
                        help="Maximum variant size to consider for comparison (-1 = off; %(default)s)")
    filteg.add_argument("--passonly", action="store_true", default=False,
                        help="Only consider calls with FILTER == PASS")

    args = parser.parse_args(args)

    return args


def check_params(args):
    """
    Checks parameters as much as possible.
    All errors are written to stderr without logging since failures mean no output
    """
    check_fail = False
    if not os.path.exists(args.input):
        check_fail = True
        logging.error("File %s does not exist", args.input)
    if not args.input.endswith(".gz"):
        check_fail = True
        logging.error("Input vcf %s does not end with .gz. Must be bgzip'd",
                      args.input)
    if not os.path.exists(args.input + '.tbi'):
        check_fail = True
        logging.error("Input vcf index %s.tbi does not exist. Must be indexed",
                      args.input)
    if args.hap and args.chain:
        check_fail = True
        logging.error("Cannot specify both --hap and --chain")
    if args.hap and args.keep != "first":
        check_fail = True
        logging.error("Using --hap must use --keep first")
    if args.reference:
        if not os.path.exists(args.reference):
            logging.error("Reference %s does not exist", args.reference)
            check_fail = True
    return check_fail


class IntraMergeOutput():
    """
    Output writer that consolidates into the first sample
    """

    def __init__(self, fn, header):
        self.fn = fn
        self.header = header
        if self.fn.endswith(".gz"):
            self.fh = gzip.GzipFile(self.fn, 'wb')
            self.isgz = True
        else:
            self.fh = open(self.fn, 'w')  # pylint: disable=consider-using-with
            self.isgz = False
        self.make_header()

    def make_header(self):
        """
        Writes intramerge header
        """
        self.header.add_line(
            '##FORMAT=<ID=SUPP,Number=1,Type=Integer,Description="Truvari collapse support flag">')
        for line in str(self.header).strip().split('\n'):
            if line.startswith("##"):
                self.write(line + '\n')
            elif line.startswith("#"):
                # new format and sample change
                data = line.strip().split('\t')[:10]
                self.write("\t".join(data) + '\n')

    def consolidate(self, entry, write_resolved=False):
        """
        Consolidate information into the first sample
        Returns a string
        """
        flag = 0
        found_gt = False
        needs_fill = set()
        for k, v in entry.samples[0].items():
            if fmt_none(v) is None:
                needs_fill.add(k)
        for idx, sample in enumerate(entry.samples):
            fmt = entry.samples[sample]
            if 1 in fmt['GT']:
                if not found_gt:
                    phased = entry.samples[idx].phased
                    entry.samples[0]['GT'] = fmt['GT']
                    entry.samples[0].phased = phased
                    found_gt = True
                flag |= 2 ** idx
            was_filled = set()
            for k in needs_fill:
                if fmt_none(fmt[k]) is not None:
                    entry.samples[0][k] = fmt[k]
                    was_filled.add(k)
            needs_fill -= was_filled
        entry.translate(self.header)
        entry.samples[0]['SUPP'] = flag
        if write_resolved:
            n_entry = entry.get_record()
            n_entry.ref = entry.get_ref()
            n_entry.alts = (entry.get_alt(),)
            entry = n_entry
        return "\t".join(str(entry).split('\t')[:10]) + '\n'

    def write(self, entry):
        """
        Writes header (str) or entries (VariantRecords)
        """
        if isinstance(entry, truvari.VariantRecord):
            entry = self.consolidate(entry, entry.params.write_resolved)
        if self.isgz:
            entry = entry.encode()
        self.fh.write(entry)

    def close(self):
        """
        Close the file handle
        """
        self.fh.close()


class CollapseOutput(dict):
    """
    Output writer for collapse
    """

    def __init__(self, args, params):
        """
        Makes all of the output files for collapse
        """
        super().__init__()
        logging.info("Params:\n%s", json.dumps(vars(args), indent=4))

        in_vcf = truvari.VariantFile(args.input, params=params)
        self["o_header"] = edit_header(in_vcf, args.median_info)
        self["c_header"] = trubench.edit_header(in_vcf)
        num_samps = len(self["o_header"].samples)
        if args.hap and num_samps != 1:
            logging.error(
                "--hap mode requires exactly one sample. Found %d", num_samps)
            sys.exit(100)
        if args.intra:
            self["output_vcf"] = IntraMergeOutput(
                args.output, self["o_header"])
        else:
            self["output_vcf"] = truvari.VariantFile(args.output, 'w',
                                                     header=self["o_header"], params=params)
        self["collap_vcf"] = truvari.VariantFile(args.removed_output, 'w',
                                                 header=self["c_header"], params=params)
        self["stats_box"] = {"collap_cnt": 0, "kept_cnt": 0,
                             "out_cnt": 0, "consol_cnt": 0}

    def write(self, collap, median_info=False):
        """
        Annotate and write kept/collapsed calls to appropriate files
        colap is a CollapsedCalls
        """
        self["stats_box"]["out_cnt"] += 1
        # Nothing collapsed, no need to annotate
        if not collap.matches:
            self["output_vcf"].write(collap.entry)
            return

        collap.annotate_entry(self["o_header"], median_info)

        self["output_vcf"].write(collap.entry)
        self["stats_box"]["kept_cnt"] += 1
        self["stats_box"]["consol_cnt"] += collap.gt_consolidate_count
        for match in collap.matches:
            trubench.annotate_entry(match.comp, match, self["c_header"])
            self["collap_vcf"].write(match.comp)
            self['stats_box']["collap_cnt"] += 1

    def close(self):
        """
        Close all the files
        """
        self["output_vcf"].close()
        self["collap_vcf"].close()

    def dump_log(self):
        """
        Log information collected during collapse
        """
        logging.info("Wrote %d Variants", self["stats_box"]["out_cnt"])
        logging.info("%d variants collapsed into %d variants",
                     self["stats_box"]["collap_cnt"], self["stats_box"]["kept_cnt"])
        logging.info("%d samples' FORMAT fields consolidated",
                     self["stats_box"]["consol_cnt"])


class LinkedList:
    """
    Simple linked list which should(?) be faster than concatenating a bunch
    regular lists
    """

    def __init__(self, data=None):
        """
        init
        """
        self.head = None
        self.tail = None
        if data is not None:
            self.append(data)

    def append(self, data):
        """
        Put data onto end of list
        """
        new_node = [data, None]
        if not self.head:
            self.head = new_node
            self.tail = new_node
            return
        self.tail[1] = new_node
        self.tail = new_node

    def to_list(self):
        """
        Turn into a regular list
        """
        cur_node = self.head
        ret = []
        while cur_node:
            ret.append(cur_node[0])
            cur_node = cur_node[1]
        return ret

    def concatenate(self, other):
        """
        Combine two linked lists
        """
        if not self.head:  # If the first list is empty
            return other
        self.tail[1] = other.head
        self.tail = other.tail
        return self


def merge_intervals(intervals):
    """
    Merge list of tuples
    """
    if not intervals:
        return []
    intervals.sort(key=lambda x: (x[0], x[1]))
    merged = []
    current_start, current_end, current_data = intervals[0]
    for i in range(1, len(intervals)):
        next_start, next_end, next_data = intervals[i]

        # Should be <=, maybe. but this replicates intervaltree
        if next_start < current_end:
            current_end = max(current_end, next_end)
            current_data.concatenate(next_data)
        else:
            # No overlap
            merged.append((current_start, current_end, current_data))
            current_start, current_end, current_data = next_start, next_end, next_data
    merged.append((current_start, current_end, current_data))
    return merged


def tree_size_chunker(params, chunks):
    """
    To reduce the number of variants in a chunk try to sub-chunk by size-similarity before hand
    Needs to return the same thing as a chunker
    """
    chunk_count = 0
    thresh = 1 if "COLLAP_SUB" in os.environ and os.environ["COLLAP_SUB"] == "1" else 100
    for chunk, _ in chunks:
        if len(chunk['base']) < thresh:  # fewer than 100 is fine
            chunk_count += 1
            yield chunk, chunk_count
            continue
        yield {'__filtered': chunk['__filtered'], 'base': []}, chunk_count
        to_add = []
        for entry in chunk['base']:
            # How much smaller/larger would be in sizesim?
            sz = entry.var_size()
            diff = sz * (1 - params.pctsize)
            if not params.typeignore:
                sz *= -1 if entry.var_type() == truvari.SV.DEL else 1
            to_add.append((sz - diff, sz + diff, LinkedList(entry)))
        tree = merge_intervals(to_add)
        for intv in tree:
            chunk_count += 1
            yield {'base': intv[2].to_list(), '__filtered': []}, chunk_count


def tree_dist_chunker(params, chunks):
    """
    To reduce the number of variants in a chunk try to sub-chunk by reference distance before hand
    Needs to return the same thing as a chunker
    This does nothing
    """
    chunk_count = 0
    thresh = 1 if "COLLAP_SUB" in os.environ and os.environ["COLLAP_SUB"] == "1" else 100
    for chunk, _ in chunks:
        if len(chunk['base']) < thresh:  # fewer than 100 is fine
            chunk_count += 1
            yield chunk, chunk_count
            continue
        yield {'__filtered': chunk['__filtered'], 'base': []}, chunk_count
        to_add = []
        for entry in chunk['base']:
            st, ed = entry.boundaries()
            st -= params.refdist
            ed += params.refdist
            to_add.append((st, ed, LinkedList(entry)))
        tree = merge_intervals(to_add)
        for intv in tree:
            chunk_count += 1
            yield {'base': intv[2].to_list(), '__filtered': []}, chunk_count


def collapse_main(args):
    """
    Main
    """
    args = parse_args(args)
    truvari.setup_logging(args.debug, show_version=True)

    if check_params(args):
        sys.stderr.write("Couldn't run Truvari. Please fix parameters\n")
        sys.exit(100)

    params = truvari.VariantParams(args=args,
                              short_circuit=True,
                              skip_gt=True,
                              sizefilt=args.sizemin,
                              chunksize=args.refdist)
    # Extra attributes needed for collapse
    params.keep = args.keep
    params.hap = args.hap
    params.gt = args.gt
    params.chain = args.chain
    params.sorter = SORTS[args.keep]
    params.no_consolidate = args.no_consolidate

    base = truvari.VariantFile(args.input, params=params)
    regions = truvari.build_region_tree(base, includebed=args.bed)
    truvari.merge_region_tree_overlaps(regions)
    base_i = base.fetch_regions(regions)

    chunks = truvari.chunker(params, ('base', base_i))
    smaller_chunks = tree_size_chunker(params, chunks)
    even_smaller_chunks = tree_dist_chunker(params, smaller_chunks)

    outputs = CollapseOutput(args, params)
    m_collap = partial(collapse_chunk, params=params)
    for call in itertools.chain.from_iterable(map(m_collap, even_smaller_chunks)):
        outputs.write(call, args.median_info)

    outputs.close()
    outputs.dump_log()
    logging.info("Finished collapse")
