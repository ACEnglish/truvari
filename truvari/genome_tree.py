"""
Helper class to specify included regions of the genome when iterating events.
"""
import re
import sys
import logging
from collections import defaultdict

from intervaltree import IntervalTree
import truvari.comparisons as tcomp

HEADERMAT=re.compile(r"##\w+=<ID=(?P<name>\w+),Number=(?P<num>[\.01AGR]),Type=(?P<type>\w+)")
class GenomeTree():
    """
    Helper class to specify included regions of the genome when iterating events.
    """
    def __init__(self, vcfA, vcfB, includebed=None, max_span=None):

        self.includebed = includebed
        self.max_span = max_span
        self.tree = self.__build_tree(vcfA, vcfB)

    def __build_tree(self, vcfA, vcfB):
        """
        Build the include regions
        """
        contigA_set = set(vcfA.header.contigs.keys())
        contigB_set = set(vcfB.header.contigs.keys())
        all_regions = defaultdict(IntervalTree)
        if self.includebed is not None:
            counter = 0
            with open(self.includebed, 'r') as fh:
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
        Iterates a vcf and yields only the entries that overlap an 'include' region
        """
        for chrom in self.tree:
            for intv in self.tree[chrom]:
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
        overlaps = self.tree[entry.chrom].overlaps(astart) and self.tree[entry.chrom].overlaps(aend)
        if astart == aend:
            return overlaps
        return overlaps and len(self.tree[entry.chrom].overlap(astart, aend)) == 1


def make_interval_tree(vcf_file, sizemin=10, sizemax=100000, passonly=False):
    """
    Return a dictionary of IntervalTree for intersection querying
    Return how many entries there are total in the vcf_file
    Return how many entries pass filtering parameters in vcf_files
    """
    n_entries = 0
    cmp_entries = 0
    lookup = defaultdict(IntervalTree)
    try:
        for entry in vcf_file:
            n_entries += 1
            if passonly and "PASS" not in entry.filter:
                continue
            start, end = tcomp.entry_boundaries(entry)
            sz = tcomp.entry_size(entry)
            if sz < sizemin or sz > sizemax:
                continue
            cmp_entries += 1
            lookup[entry.chrom].addi(start, end, entry.start)
    except ValueError as e:
        logging.error("Unable to parse comparison vcf file. Please check header definitions")
        logging.error("Specific error: \"%s\"", str(e))
        sys.exit(100)

    return lookup, n_entries, cmp_entries
