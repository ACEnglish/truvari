"""
Helper class to specify included regions of the genome when iterating events.
"""
import sys
import copy
import logging
from collections import defaultdict

from intervaltree import IntervalTree
import truvari
import truvari.comparisons as tcomp


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
                "Excluding %d contigs present in comparison calls header but not base calls.", len(excluding))

        for contig in contigA_set:
            name = vcfA.header.contigs[contig].name
            length = vcfA.header.contigs[contig].length
            if not length:
                logging.error("Contig %s has no length definition. Fix header.", name)
                sys.exit(10)
            all_regions[name].addi(0, length)
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
                    logging.warning("Unable to fetch %s from %s", chrom, vcf_file.filename)

    def include(self, entry):
        """
        Returns if this entry's start and end are within a region that is to be included
        Here overlap means lies completely within the boundary of an include region
        """
        astart, aend = tcomp.entry_boundaries(entry)
        # Filter these early so we don't have to keep checking overlaps
        if self.max_span is None or aend - astart > self.max_span:
            return False
        overlaps = self.tree[entry.chrom].overlaps(astart) \
            and self.tree[entry.chrom].overlaps(aend)
        if astart == aend:
            return overlaps
        return overlaps and len(self.tree[entry.chrom].overlap(astart, aend)) == 1

    def extend(self, pad):
        """
        Extends all intervals by a fixed number of bases on each side
        Returns a copy of this IntervalTree
        """
        logging.info("Extending the regions by %d bases", pad)
        ret = copy.deepcopy(self)

        for chrom in ret.tree:
            ret.tree[chrom] = IntervalTree.from_tuples(((max(0, i.begin - pad), i.end + pad)) for i in ret.tree[chrom])

        ret.merge_overlaps()
        return ret

def build_anno_tree(filename, chrom_col=0, start_col=1, end_col=2, one_based=False, comment='#', idxfmt=None):
    """
    Build an dictionary of IntervalTrees for each chromosome from tab-delimited annotation file
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
        tree[chrom].addi(start, end, data=m_idx)
        idx += 1
    return tree, idx
