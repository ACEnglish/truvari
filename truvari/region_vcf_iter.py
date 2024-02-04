"""
Helper class to specify included regions of the genome when iterating events.
"""
import sys
import copy
import logging
from collections import defaultdict, deque

from intervaltree import IntervalTree
import truvari

class RegionVCFIterator():
    """
    Helper class to specify include regions of the genome when iterating a VCF
    Subset to only events less than max_span.
    Subset to only events on contigs listed in vcfA left-join vcfB
    """

    def __init__(self, vcfA, vcfB=None, includebed=None, max_span=None):
        """ init """
        self.includebed = includebed
        self.max_span = max_span
        self.tree = self.__build_tree(vcfA, vcfB)

    def __build_tree(self, vcfA, vcfB):
        """
        Build the include regions
        """
        contigA_set = set(vcfA.header.contigs.keys())
        if vcfB is not None:
            contigB_set = set(vcfB.header.contigs.keys())
        else:
            contigB_set = contigA_set

        if self.includebed is not None:
            all_regions, counter = build_anno_tree(self.includebed)
            logging.info("Including %d bed regions", counter)
            return all_regions

        all_regions = defaultdict(IntervalTree)
        excluding = contigB_set - contigA_set
        if excluding:
            logging.warning(
                "Excluding %d contigs present in comparison calls header but not baseline calls.", len(excluding))

        for contig in contigA_set:
            name = vcfA.header.contigs[contig].name
            length = vcfA.header.contigs[contig].length
            if not length:
                logging.error(
                    "Contig %s has no length definition. Fix header.", name)
                sys.exit(10)
            all_regions[name].addi(0, length + 1)
        return all_regions

    def merge_overlaps(self):
        """
        Runs IntervalTree.merge_overlaps on all trees. Returns list of all chromosomes having overlapping regions
        and the pre_post totals
        """
        chr_with_overlaps = []
        for i in self.tree:
            pre_len = len(self.tree[i])
            self.tree[i].merge_overlaps()
            post_len = len(self.tree[i])
            if pre_len != post_len:
                chr_with_overlaps.append(i)
        if chr_with_overlaps:
            logging.info("Found %d chromosomes with overlapping regions",
                         len(chr_with_overlaps))
            logging.debug("CHRs: %s", chr_with_overlaps)

    def iterate(self, vcf_file):
        """
        Iterates a vcf and yields only the entries that overlap included regions
        """
        for chrom in sorted(self.tree.keys()):
            for intv in sorted(self.tree[chrom]):
                try:
                    for entry in vcf_file.fetch(chrom, intv.begin, intv.end):
                        if self.includebed is None or self.include(entry):
                            yield entry
                except ValueError:
                    logging.warning("Unable to fetch %s from %s",
                                    chrom, vcf_file.filename)

    def include(self, entry):
        """
        Returns if this entry's start and end are within a region that is to be included
        Here overlap means lies completely within the boundary of an include region
        """
        astart, aend = truvari.entry_boundaries(entry)
        # Filter these early so we don't have to keep checking overlaps
        if self.max_span is not None and aend - astart > self.max_span:
            return False

        m_ovl = self.tree[entry.chrom].overlap(astart, aend)
        if len(m_ovl) != 1:
            return False
        m_ovl = list(m_ovl)[0]
        end_within = truvari.entry_variant_type(entry) != truvari.SV.INS
        return truvari.coords_within(astart, aend, m_ovl.begin, m_ovl.end - 1, end_within)

    def extend(self, pad):
        """
        Extends all intervals by a fixed number of bases on each side
        Returns a copy of this IntervalTree
        """
        logging.info("Extending the regions by %d bases", pad)
        ret = copy.deepcopy(self)

        for chrom in ret.tree:
            ret.tree[chrom] = IntervalTree.from_tuples(
                ((max(0, i.begin - pad), i.end + pad)) for i in ret.tree[chrom])

        ret.merge_overlaps()
        return ret


def build_anno_tree(filename, chrom_col=0, start_col=1, end_col=2, one_based=False, comment='#', idxfmt=None):
    """
    Build an dictionary of IntervalTrees for each chromosome from tab-delimited annotation file

    By default, the file is assumed to be a bed-format. If custom chrom/start/end are used, the columns can be
    specified.

    idxfmt is a string that will be formatted with the line number from the file (e.g. "num %d" will
    make intervals with data="num 0", "num 1". By default the data will be interger line number.
    If intervals will be compared between anno_trees, set idxfmt to ""

    :param `filename`: Path to file to parse, can be compressed
    :type `filename`: string
    :param `chrom_col`: Index of column in file with chromosome
    :type `chrom_col`: integer, optional
    :param `start_col`: Index of column in file with start
    :type `start_col`: integer, optional
    :param `end_col`: Index of column in file with end
    :type `end_col`: integer, optional
    :param `one_based`: True if coordinates are one-based
    :type `one_based`: bool, optional
    :param `comment`: ignore lines if they start with this string
    :type `comment`: string, optional
    :param `idxfmt`: Index of column in file with chromosome
    :type `idxfmt`: string, optional

    :return: dictionary with chromosome keys and :class:`intervaltree.IntervalTree` values
    :rtype: dict
    """
    idx = 0
    correction = 1 if one_based else 0
    tree = defaultdict(IntervalTree)
    for line in truvari.opt_gz_open(filename):
        if line.startswith(comment):
            continue
        data = line.strip().split('\t')
        chrom = data[chrom_col]
        start = int(data[start_col]) - correction
        end = int(data[end_col])
        if idxfmt is not None:
            m_idx = idxfmt.format(idx)
        else:
            m_idx = idx
        tree[chrom].addi(start, end + 1, data=m_idx)
        idx += 1
    return tree, idx

def region_filter(vcf, tree, inside=True):
    """
    Given a VariantRecord iter and defaultdict(IntervalTree), yield variants which are inside/outside the tree regions
    """
    for chrom, cur_tree in tree.items():
        cur_tree = deque(sorted(cur_tree))
        try:
            cur_intv = cur_tree.popleft()
        except IndexError as e:
            # region-less chromosome
            continue

        for entry in vcf.fetch(chrom):
            start, end = truvari.entry_boundaries(entry)
            # if start is after this interval, we need the next interval
            if start > cur_intv.end - 1:
                try:
                    cur_intv = cur_tree.popleft()
                except IndexError as e:
                    if not inside:
                        yield entry
                        continue
                    # next chromosome
                    break 

            # well before
            if end < cur_intv.begin:
                if not inside:
                    yield entry
                continue

            end_within = truvari.entry_variant_type(entry) != truvari.SV.INS
            is_within = truvari.coords_within(start, end, cur_intv.begin, cur_intv.end - 1, end_within)
            if is_within == inside:
                yield entry

