"""
Helper class to specify included regions of the genome when iterating events.
"""
import io
import re
import sys
import gzip
import logging
from collections import defaultdict, OrderedDict

from intervaltree import IntervalTree

import truvari.comparisons as tcomp

HEADERMAT=re.compile("##\w+=<ID=(?P<name>\w+),Number=(?P<num>[\.01AGR]),Type=(?P<type>\w+)")
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
            if not excluding:
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

def make_bedanno_tree(bed_file):
    """
    Expecting a bedfile formatted with new vcf headerlines
    at the top then bed entries the rest of the way
    Okay = you need VCF HEADER lines at the top
    For now, this will only annotate INFO fields
    Those header lines must be in order of the bed columns' annotations
    So as an example:

    ##INFO=<ID=SREP,Number=0,Type=Flag,Description="Simple Repeat">
    ##INFO=<ID=SREPLEN,Number=1,Type=Integer,Description="Simple Repeat Length">
    1	10000	10468	trf	789
    1	10627	10800	.	346
    
    This will add the SREP flag IF column[3] != '.'
    This will add the SREPLEN as an Integer from column[4]
    Errors will be thrown if there's a problem parsing

    returns a tuple of the (defaultdict(IntervalTree), header_lines)
    """
    logging.info("Loading Annotation %s", bed_file)
    n_entries = 0
    lookup = defaultdict(IntervalTree)
    header_lines = []
    header_dict = OrderedDict()

    def add_header(line):
        header_lines.append(line)
        head = HEADERMAT.match(line)
        if not head:
            logging.error("Bad Headerline %s", line.strip())
            exit(1)
        g = head.groupdict()
        typ = None
        if g["type"] == "String":
            typ = str
        elif g["type"] == "Integer":
            typ = int
        elif g["type"] == "Float":
            typ = float
        elif g["type"] == "Flag":
            typ = bool
        
        if not typ:
            logging.error("Bad Headerline Type %s", line.strip())
            exit(1)
                
        num = None
        if g["num"] == '0':
            num = bool
        elif g["num"] != 1:
            num = list
        else:
            num = typ
        header_dict[g["name"]] = (num, typ)

    # I think I'm going to need to build a header parser that will convert the things
    if bed_file.endswith(".gz"):
        fh = io.TextIOWrapper(gzip.open(bed_file))
    else:
        fh = open(bed_file, 'r')
    
    for line in fh:
        # header lines - ##\w+=<ID=(?P<name>\w+),Number=(?P<num>[\.01AGR]),Type=(?P<type>\w+)
        if line.startswith("#"):
            add_header(line)
            continue
        
        # parse the rest
        data = line.strip().split('\t')
        start = int(data[1])
        end = int(data[2])
        if len(data[3:]) != len(header_dict):
            logging.error("Bad header in file %s. %d headers and %d columns on line %s", 
                          bed_file, len(header_dict), len(data[3:]), line)
            exit(1)
        m_dict = {}
        for k, v in zip(header_dict.keys(), data[3:]):
            if v != ".":
                if header_dict[k][0] == list:
                    m_dict[k] = [header_dict[k][1](x) for x in v.split(',')]
                elif header_dict[k][0] == bool: 
                    m_dict[k] = True
                else:
                    m_dict[k] = header_dict[k][1](v)
        lookup[data[0]].addi(start, end, m_dict)
    return lookup, header_lines
