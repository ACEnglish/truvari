"""
Calcluate the the Hom / (Het + Hom) of variants in the region of SVs
Requires the VCF to contain SVs beside SNPs/Indels
"""

import sys
import logging
import argparse

import pysam
import truvari
from acebinf import setup_logging

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="hompct", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Compressed, indexed VCF to annotate ")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output filename (stdout)")
    parser.add_argument("-b", "--buffer", type=int, default=5000,
                        help="Number of base-pairs up/dn-stream to query (%(default)s)")
    parser.add_argument("-m", "--minanno", type=int, default=50,
                        help="Minimum size of event to annotate")
    parser.add_argument("-M", "--maxgt", type=int, default=1,
                        help="Largest event size to count for genotyping (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    setup_logging(args.debug)
    return args

def hompct_main(cmd_args):
    args = parse_args(cmd_args)
    
    v = pysam.VariantFile(args.input)
    def get_pct(chrom, start, end):
        tot = 0
        homs = 0
        for entry in v.fetch(chrom, max(0, start - args.buffer), min(v.header.contigs[chrom].length, end + args.buffer)):
            if truvari.entry_size(entry) > args.maxgt:
                continue
            if truvari.get_gt(entry.samples[0]["GT"]).name == "HOM":
                homs += 1
            tot += 1
                
        if tot == 0:
            return float('nan')

        return float(format((homs / tot) * 100, ".1f"))
            
    header = v.header.copy()
    header.add_line(('##INFO=<ID=HOMPCT,Number=1,Type=Float,'
                     'Description="Percent of calls < %dbp long within %dbp that are homozygous') 
                     % (args.maxgt, args.buffer))

    out = pysam.VariantFile(args.output, 'w', header=header)
    v2 = pysam.VariantFile(args.input)
    for entry in v2:
        if truvari.entry_size(entry) >= args.minanno:
            entry = truvari.copy_entry(entry, header)
            entry.info["HOMPCT"] = get_pct(entry.chrom, *truvari.entry_boundaries(entry))
        out.write(entry)
    logging.info("Finished hompct")

if __name__ == '__main__':
    hompct_main(sys.argv[1:])
