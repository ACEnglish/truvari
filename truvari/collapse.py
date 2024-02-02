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

import pysam
import numpy as np
from intervaltree import IntervalTree

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
    genotype_mask: str = ""  # bad

    def combine(self, other):
        """
        Put other's entries into this' collapsed_entries
        """
        self.matches.append(other.entry)
        self.matches.extend(other.matches)

    def calc_median_sizepos(self):
        """
        Given a set of variants, calculate the median start, end, and size
        """
        st, en = truvari.entry_boundaries(self.entry)
        sz = truvari.entry_size(self.entry)
        starts = [st]
        ends = [en]
        sizes = [sz]
        for i in self.matches:
            st, en = truvari.entry_boundaries(i.comp)
            sz = truvari.entry_size(i.comp)
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


def chain_collapse(cur_collapse, all_collapse, matcher):
    """
    Perform transitive matching of cur_collapse to all_collapse
    """
    for m_collap in all_collapse:
        for other in m_collap.matches:
            mat = matcher.build_match(cur_collapse.entry,
                                      other.base,
                                      m_collap.match_id,
                                      skip_gt=True,
                                      short_circuit=True)
            if mat.state:
                m_collap.consolidate(cur_collapse)
                return True  # you can just ignore it later
    return False  # You'll have to add it to your set of collapsed calls


def collapse_chunk(chunk, matcher):
    """
    Returns a list of lists with [keep entry, collap matches, match_id, gt_consolidate_count]
    """
    chunk_dict, chunk_id = chunk
    remaining_calls = sorted(chunk_dict['base'], key=matcher.sorter)

    remaining_calls.sort(key=matcher.sorter)
    call_id = -1
    ret = []  # list of Collapses
    while remaining_calls:
        call_id += 1
        m_collap = CollapsedCalls(remaining_calls.pop(0),
                                  f'{chunk_id}.{call_id}')
        # quicker genotype comparison - needs to be refactored
        m_collap.genotype_mask = m_collap.make_genotype_mask(
            m_collap.entry, matcher.gt)

        # Sort based on size difference to current call
        for candidate in sorted(remaining_calls, key=partial(relative_size_sorter, m_collap.entry)):
            mat = matcher.build_match(m_collap.entry,
                                      candidate,
                                      m_collap.match_id,
                                      skip_gt=True,
                                      short_circuit=True)
            if matcher.hap and not hap_resolve(m_collap.entry,  candidate):
                mat.state = False
            if mat.state and m_collap.gt_conflict(candidate, matcher.gt):
                mat.state = False
            if mat.state:
                m_collap.matches.append(mat)

        # Does this collap need to go into a previous collap?
        if not matcher.chain or not chain_collapse(m_collap, ret, matcher):
            ret.append(m_collap)

        # If hap, only allow the best match
        if matcher.hap and m_collap.matches:
            mats = sorted(m_collap.matches, reverse=True)
            m_collap.matches = [mats.pop(0)]

        # Remove everything that was used
        to_rm = [_.comp for _ in m_collap.matches]
        remaining_calls = [_ for _ in remaining_calls if _ not in to_rm]

    if matcher.no_consolidate:
        for val in ret:
            if matcher.gt != 'off':
                edited_entry, collapse_cnt = gt_aware_consolidate(
                    val.entry, val.matches)
            else:
                edited_entry, collapse_cnt = collapse_into_entry(
                    val.entry, val.matches, matcher.hap)
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
    return abs(truvari.entry_size(base) - truvari.entry_size(comp))


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
                except TypeError:
                    # Happens for things like PL when one is null but its expecting a tuple
                    logging.debug("Unable to set FORMAT %s for sample %s",
                                  key, sample)
                    logging.debug("Kept entry: %s:%d %s",
                                  entry.chrom, entry.pos, entry.id)
                    logging.debug("Colap entry: %s:%d %s",
                                  o_entry.chrom, o_entry.pos, o_entry.id)
                except KeyError:
                    logging.debug("Unshared format %s in sample %s ignored for pair %s:%d %s %s:%d %s",
                                  key, sample, entry.chrom, entry.pos, entry.id, o_entry.chrom,
                                  o_entry.pos, o_entry.id)
            # pass along phase
            entry.samples[sample].phased = o_entry.samples[sample].phased
    return entry, n_consolidate


def gt_conflict(cur_collapse, entry, which_gt):
    """
    Return true if entry's genotypes conflict with any of the current collapse
    which_gt all prevents variants present in the same sample from being collapsed
    which_gt het only prevents two het variants from being collapsed.
    Might be deprecated, now?
    """
    if which_gt == 'off':
        return False

    def checker(base, comp):
        """
        Return true if a conflict
        """
        i = 0
        samps = list(base.samples)
        while i < len(samps):
            sample = samps[i]
            gtA = base.samples[sample].allele_indices
            gtB = comp.samples[sample].allele_indices
            if which_gt == 'all':
                if 1 in gtA and 1 in gtB:
                    return True
            elif gtA.count(1) == 1 and gtB.count(1) == 1:
                return True
            i += 1
        return False
    # Need to check the kept entry
    if checker(cur_collapse.entry, entry):
        return True
    # And all removed entries
    for mat in cur_collapse.matches:
        if checker(entry, mat.comp):
            return True

    return False


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
    if isinstance(value, tuple):
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
            new_fmts[k] = None if k not in entry.samples.keys(
            ) else entry.samples[sample][k]
        for o in others:
            o = o.comp
            n_c = get_ac(o.samples[sample]['GT'])
            n_consolidated += 1 if n_c != 0 else 0
            new_fmts['GT'] += n_c
            for k in all_fmts:
                rpl_val = None if k not in o.samples[sample] else fmt_none(
                    o.samples[sample][k])
                if not (k in new_fmts and new_fmts[k] is not None) and rpl_val:
                    new_fmts[k] = o.samples[sample][k]
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
    s1 = truvari.entry_size(b1)
    s2 = truvari.entry_size(b2)
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
    Order entries from highest to lowest MAC
    """
    mac1 = truvari.allele_freq_annos(b1)["MAC"]
    mac2 = truvari.allele_freq_annos(b2)["MAC"]
    if mac1 < mac2:
        return 1
    if mac1 > mac2:
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
    parser.add_argument("-c", "--collapsed-output", type=str, default="collapsed.vcf",
                        help="Where collapsed variants are written (collapsed.vcf)")
    parser.add_argument("-f", "--reference", type=str, default=None,
                        help="Indexed fasta used to call variants")
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
    thresg.add_argument("-B", "--minhaplen", type=truvari.restricted_int, default=50,
                        help="Minimum haplotype sequence length to create (%(default)s)")
    thresg.add_argument("-P", "--pctsize", type=truvari.restricted_float, default=0.95,
                        help="Min pct allele size similarity (minvarsize/maxvarsize) (%(default)s)")
    thresg.add_argument("-O", "--pctovl", type=truvari.restricted_float, default=0.0,
                        help="Min pct reciprocal overlap (%(default)s) for DEL events")
    thresg.add_argument("-t", "--typeignore", action="store_true", default=False,
                        help="Variant types don't need to match to compare (%(default)s)")

    parser.add_argument("--hap", action="store_true", default=False,
                        help="Collapsing a single individual's haplotype resolved calls (%(default)s)")
    parser.add_argument("--chain", action="store_true", default=False,
                        help="Chain comparisons to extend possible collapsing (%(default)s)")
    parser.add_argument("--no-consolidate", action="store_false", default=True,
                        help="Skip consolidation of sample genotype fields (%(default)s)")
    parser.add_argument("--null-consolidate", type=str, default=None,
                        help=("Comma separated list of FORMAT fields to consolidate into the kept "
                              "entry by taking the first non-null from all neighbors (%(default)s)"))
    filteg = parser.add_argument_group("Filtering Arguments")
    filteg.add_argument("-s", "--sizemin", type=truvari.restricted_int, default=50,
                        help="Minimum variant size to consider for comparison (%(default)s)")
    filteg.add_argument("-S", "--sizemax", type=truvari.restricted_int, default=50000,
                        help="Maximum variant size to consider for comparison (%(default)s)")
    filteg.add_argument("--passonly", action="store_true", default=False,
                        help="Only consider calls with FILTER == PASS")

    args = parser.parse_args(args)

    if args.null_consolidate is not None:
        args.null_consolidate = args.null_consolidate.split(',')

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

    def consolidate(self, entry):
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
        return "\t".join(str(entry).split('\t')[:10]) + '\n'

    def write(self, entry):
        """
        Writes header (str) or entries (VariantRecords)
        """
        if isinstance(entry, pysam.VariantRecord):
            entry = self.consolidate(entry)
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

    def __init__(self, args):
        """
        Makes all of the output files for collapse
        """
        super().__init__()
        logging.info("Params:\n%s", json.dumps(vars(args), indent=4))

        in_vcf = pysam.VariantFile(args.input)
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
            self["output_vcf"] = pysam.VariantFile(args.output, 'w',
                                                   header=self["o_header"])
        self["collap_vcf"] = pysam.VariantFile(args.collapsed_output, 'w',
                                               header=self["c_header"])
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


def tree_size_chunker(matcher, chunks):
    """
    To reduce the number of variants in a chunk try to sub-chunk by size-similarity before hand
    Needs to return the same thing as a chunker
    """
    chunk_count = 0
    for chunk, _ in chunks:
        if len(chunk['base']) < 100:  # fewer than 100 is fine
            chunk_count += 1
            yield chunk, chunk_count
            continue
        tree = IntervalTree()
        yield {'__filtered': chunk['__filtered'], 'base': []}, chunk_count
        for entry in chunk['base']:
            # How much smaller/larger would be in sizesim?
            sz = truvari.entry_size(entry)
            diff = sz * (1 - matcher.params.pctsize)
            if not matcher.params.typeignore:
                sz *= - \
                    1 if truvari.entry_variant_type(
                        entry) == truvari.SV.DEL else 1
            tree.addi(sz - diff, sz + diff, data=[entry])
        tree.merge_overlaps(data_reducer=lambda x, y: x + y)
        for intv in tree:
            chunk_count += 1
            yield {'base': intv.data, '__filtered': []}, chunk_count


def tree_dist_chunker(matcher, chunks):
    """
    To reduce the number of variants in a chunk try to sub-chunk by reference distance before hand
    Needs to return the same thing as a chunker
    """
    chunk_count = 0
    for chunk, _ in chunks:
        if len(chunk['base']) < 100:  # fewer than 100 is fine
            chunk_count += 1
            yield chunk, chunk_count
            continue
        tree = IntervalTree()
        yield {'__filtered': chunk['__filtered'], 'base': []}, chunk_count
        for entry in chunk['base']:
            st, ed = truvari.entry_boundaries(entry)
            st -= matcher.params.refdist
            ed += matcher.params.refdist
            tree.addi(st, ed, data=[entry])
        tree.merge_overlaps(data_reducer=lambda x, y: x + y)
        for intv in tree:
            chunk_count += 1
            yield {'base': intv.data, '__filtered': []}, chunk_count


def collapse_main(args):
    """
    Main
    """
    args = parse_args(args)
    truvari.setup_logging(args.debug, show_version=True)

    if check_params(args):
        sys.stderr.write("Couldn't run Truvari. Please fix parameters\n")
        sys.exit(100)

    # The matcher collapse needs isn't 100% identical to bench
    # So we set it up here
    args.chunksize = args.refdist
    args.bSample = None
    args.cSample = None
    args.sizefilt = args.sizemin
    args.no_ref = False
    matcher = truvari.Matcher(args=args)
    matcher.params.includebed = None
    matcher.keep = args.keep
    matcher.hap = args.hap
    matcher.gt = args.gt
    matcher.chain = args.chain
    matcher.sorter = SORTS[args.keep]
    matcher.no_consolidate = args.no_consolidate
    matcher.picker = 'single'

    base = pysam.VariantFile(args.input)
    regions = truvari.RegionVCFIterator(base, includebed=args.bed)
    regions.merge_overlaps()
    base_i = regions.iterate(base)

    chunks = truvari.chunker(matcher, ('base', base_i))
    smaller_chunks = tree_size_chunker(matcher, chunks)
    even_smaller_chunks = tree_dist_chunker(matcher, smaller_chunks)

    outputs = CollapseOutput(args)
    m_collap = partial(collapse_chunk, matcher=matcher)
    for call in itertools.chain.from_iterable(map(m_collap, even_smaller_chunks)):
        outputs.write(call, args.median_info)

    outputs.close()
    outputs.dump_log()
    logging.info("Finished collapse")
