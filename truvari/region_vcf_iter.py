"""
Helper class to specify included regions of the genome when iterating events.
"""
import logging
from collections import defaultdict

from intervaltree import IntervalTree
import truvari.comparisons as tcomp
from truvari.annos.bpovl import build_anno_tree

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
            all_regions[name].addi(0, length)
        return all_regions

    def iterate(self, vcf_file):
        """
        Iterates a vcf and yields only the entries that overlap included regions
        """
        for chrom in sorted(self.tree.keys()):
            for intv in sorted(self.tree[chrom]):
                for entry in vcf_file.fetch(chrom, intv.begin, intv.end):
                    if self.includebed is None or self.include(entry):
                        yield entry

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
