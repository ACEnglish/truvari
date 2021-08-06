"""
For every call within size boundaries,
Add NumNeighbors info field of how many calls are within the distance
Add NeighId clustering field in the same chained neighborhood
For example,
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
    parser.add_argument("-r", "--refdist", type=int, default=1000,
                        help="Max reference location distance (%(default)s)")
    parser.add_argument("-s", "--sizemin", type=int, default=50,
                        help="Minimum variant size to consider for annotation (%(default)s)")
    parser.add_argument("--passonly", action="store_true", default=False,
                        help="Only count calls with FILTER == PASS")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    return args

def edit_header(my_vcf):
    """
    Add INFO for new fields to vcf
    """
    header = my_vcf.header.copy()
    header.add_line(('##INFO=<ID=NumNeighbors,Number=1,Type=Integer,'
                     'Description="Number of calls in the neighborhood of this call">'))
    header.add_line(('##INFO=<ID=NeighId,Number=1,Type=Integer,'
                     'Description="Identifier of calls in the same chained neighborhood">'))

    return header

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

    in_vcf = pysam.VariantFile(args.input)
    header = edit_header(in_vcf)
    out_vcf = pysam.VariantFile(args.output, 'w', header=header)
    BUF = args.refdist

    def make_range(entry):
        start, end = truvari.entry_boundaries(entry)
        return [start, end, entry, 0]

    def overlaps(range1, range2):
        if range1[2].chrom != range2[2].chrom:
            return False
        start = max(0, range1[0] - BUF - 1)
        end = range1[1] + BUF + 1
        return truvari.overlaps(start, end, *range2[:2])

    NEIGHID = 0
    def output(entry, neigh_cnt):
        new_entry = truvari.copy_entry(entry, header)
        new_entry.info["NumNeighbors"] = neigh_cnt
        new_entry.info["NeighID"] = NEIGHID
        out_vcf.write(new_entry)

    def flush_push_stack(cur_range, stack):
        while stack and not overlaps(cur_range, stack[0]):
            # We know the first thing in the stack has all its downstream neighbors
            # currently in the stack and it's been updated with itps' upstream neighbors,
            # So output that one
            to_output = stack.pop(0)

            to_output[3] += len(stack)
            # the entry and it's count
            output(to_output[2], to_output[3])
            # If this event is not extending the stack, we're at the next 'locus'
            if not stack:
                NEIGHID += 1
            # Let the downstream vars know about their now flushed neighbor
            for i in stack:
                i[3] += 1
        stack.append(cur_range)

    def chrom_end_flush(stack):
        # final stack flush
        for pos, i in enumerate(stack):
            for j in stack[pos + 1:]:
                if overlaps(i, j):
                    i[3] += 1
                if overlaps(j, i):
                    j[3] += 1

        for i in stack:
            output(i[2], i[3])

    stack = []
    last_pos = None
    for entry in in_vcf:
        size = truvari.entry_size(entry)
        if not last_pos:
            last_pos = [entry.chrom, entry.start]

        if last_pos[0] == entry.chrom and last_pos[1] > entry.start:
            logging.error("File is not sorted %d:%d before %s:%d",
                          last_pos[0], last_pos[1], entry.chrom, entry.start)
            sys.exit(1)

        if entry.chrom != last_pos[0]:
            chrom_end_flush(stack)
            stack = []

        last_pos = [entry.chrom, entry.start]

        if size < args.sizemin or (args.passonly and truvari.filter_value(entry)):
            out_vcf.write(entry)
            continue

        cur_range = make_range(entry)
        flush_push_stack(cur_range, stack)

    chrom_end_flush(stack)

    logging.info("Finished")

if __name__ == '__main__':
    numneigh_main(sys.argv[1:])
