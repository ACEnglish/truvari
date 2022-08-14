"""
For every call within size boundaries,
Add NumNeighbors info field of how many calls are within the distance
Add NeighId clustering field in the same chained neighborhood
For example,
::
    -- is a call, refdist is 2
         - - -   -    - -
    nn:  1 2 1   0    1 1
    id:  0 0 0   1    2 2

"""
import os
import sys
import logging
import argparse

import pysam
import truvari


def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="numneigh", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, default="/dev/stdin",
                        help="VCF to annotate")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output vcf (stdout)")
    parser.add_argument("-r", "--refdist", type=truvari.restricted_int, default=1000,
                        help="Max reference location distance (%(default)s)")
    parser.add_argument("-s", "--sizemin", type=truvari.restricted_int, default=50,
                        help="Minimum variant size to consider for annotation (%(default)s)")
    parser.add_argument("--passonly", action="store_true", default=False,
                        help="Only count calls with FILTER == PASS")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    return args


class NeighAnno():
    """
    Annotates a vcf with Neighbor information
    """

    def __init__(self, in_vcf, out_vcf, refdist=1000, sizemin=50, passonly=False):
        """
        initialize
        """
        self.in_vcf = pysam.VariantFile(in_vcf)
        self.header = self.edit_header()
        self.out_vcf = pysam.VariantFile(out_vcf, 'w', header=self.header)
        self.refdist = refdist
        self.sizemin = sizemin
        self.passonly = passonly
        self.neigh_id = 0
        self.stack = []

    def edit_header(self):
        """
        Add INFO for new fields to vcf
        """
        header = self.in_vcf.header.copy()
        header.add_line(('##INFO=<ID=NumNeighbors,Number=1,Type=Integer,'
                         'Description="Number of calls in the neighborhood of this call">'))
        header.add_line(('##INFO=<ID=NeighId,Number=1,Type=Integer,'
                         'Description="Identifier of calls in the same chained neighborhood">'))
        return header

    def overlaps(self, range1, range2):
        """
        Check if two ranges overlap
        """
        if range1[2].chrom != range2[2].chrom:
            return False
        start = max(0, range1[0] - self.refdist - 1)
        end = range1[1] + self.refdist + 1
        return truvari.overlaps(start, end, *range2[:2])

    def output(self, entry, neigh_cnt):
        """
        Annotate and write an entry
        """
        entry.translate(self.header)
        entry.info["NumNeighbors"] = neigh_cnt
        entry.info["NeighId"] = self.neigh_id
        self.out_vcf.write(entry)

    def flush_push_stack(self, cur_range):
        """
        Pop/output entries in stack
        """
        while self.stack and not self.overlaps(cur_range, self.stack[0]):
            # We know the first thing in the stack has all its downstream neighbors
            # currently in the stack and it's been updated with itps' upstream neighbors,
            # So output that one
            to_output = self.stack.pop(0)

            to_output[3] += len(self.stack)
            # the entry and it's count
            self.output(to_output[2], to_output[3])
            # If this event is not extending the stack, we're at the next 'locus'
            if not self.stack:
                self.neigh_id += 1
            # Let the downstream vars know about their now flushed neighbor
            for i in self.stack:
                i[3] += 1
        self.stack.append(cur_range)

    def chrom_end_flush(self):
        """
        Pop/output all entries in a stack
        """
        # final stack flush
        for pos, i in enumerate(self.stack):
            for j in self.stack[pos + 1:]:
                if self.overlaps(i, j):
                    i[3] += 1
                if self.overlaps(j, i):
                    j[3] += 1

        for i in self.stack:
            self.output(i[2], i[3])
        self.stack = []

    def run(self):
        """
        Find neighbors through a vcf
        """
        last_pos = None
        for entry in self.in_vcf:
            size = truvari.entry_size(entry)
            if not last_pos:
                last_pos = [entry.chrom, entry.start]

            if last_pos[0] == entry.chrom and last_pos[1] > entry.start:
                logging.error("File is not sorted %s:%d before %s:%d",
                              last_pos[0], last_pos[1], entry.chrom, entry.start)
                sys.exit(1)

            if entry.chrom != last_pos[0]:
                self.chrom_end_flush()
                self.neigh_id += 1

            last_pos = [entry.chrom, entry.start]

            if size < self.sizemin or (self.passonly and truvari.filter_value(entry)):
                self.out_vcf.write(entry)
                continue
            # Make new range
            start, end = truvari.entry_boundaries(entry)
            cur_range = [start, end, entry, 0]

            self.flush_push_stack(cur_range)

        self.chrom_end_flush()


def numneigh_main(args):
    """
    Main
    """
    args = parse_args(args)
    check_fail = False
    if not os.path.exists(args.input):
        check_fail = True
        logging.error("File %s does not exist", args.input)
    if check_fail:
        logging.error("Error parsing input. Please fix before re-running")
        sys.exit(100)

    anno = NeighAnno(args.input, args.output, args.refdist,
                     args.sizemin, args.passonly)
    anno.run()
    logging.info("Finished numneigh")
