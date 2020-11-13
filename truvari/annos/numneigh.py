"""
For every call within size boundaries, 
Add NumNeighbors info field of how many calls over within size boundary
are in the neighborhood
"""
import os
import logging
import argparse

import pysam
import truvari

def parse_args(args):
    """
    Pull the command line parameters
    """
    def restricted_float(x):
        x = float(x)
        if x < 0.0 or x > 1.0:
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
        return x

    parser = argparse.ArgumentParser(prog="numneigh", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, required=True,
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
    return header

def numneigh_main(args):
    """
    Main
    """
    args = parse_args(args)
    # Make sure we're working with an existing, compressed, indexed file
    check_fail = False
    if not os.path.exists(args.input):
        check_fail = True
        logging.error("File %s does not exist", args.input)
    if not args.input.endswith(".gz"):
        check_fail = True
        logging.error("Input vcf %s does not end with .gz. Must be bgzip'd", args.input)
    if not os.path.exists(args.input + '.tbi'):
        check_fail = True
        logging.error("Input vcf index %s.tbi does not exist. Must be indexed", args.input)
    if check_fail:
        logging.error("Error parsing input. Please fix before re-running")
        exit(100)

    in_vcf = pysam.VariantFile(args.input)
    seek_vcf = pysam.VariantFile(args.input)
    header = edit_header(in_vcf)
    out_vcf = pysam.VariantFile(args.output, 'w', header=header)

    for entry in in_vcf:
        neigh_cnt = 0
        size = truvari.entry_size(entry)
        if size < args.sizemin or (args.passonly and truvari.filter_value(entry)):
            out_vcf.write(entry)
            continue
        start, end = truvari.entry_boundaries(entry)
        start = max(0, start - args.refdist - 1)
        end = end + args.refdist + 1
        for neigh in seek_vcf.fetch(entry.chrom, start, end):
            size = truvari.entry_size(neigh)
            if not (size < args.sizemin or (args.passonly and truvari.filter_value(neigh))):
                neigh_cnt += 1
        
        # Call can't be a neighbor of itself
        neigh_cnt -= 1
        new_entry = truvari.copy_entry(entry, header)
        new_entry.info["NumNeighbors"] = neigh_cnt
        out_vcf.write(new_entry)
    logging.info("Finished")
