"""
Adds SVTYPE and SVLEN INFO fields
"""
import logging
import argparse

import pysam
import truvari


def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="svinfo", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, default="/dev/stdin",
                        help="VCF to annotate (stdin)")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output filename (stdout)")
    parser.add_argument("-m", "--minsize", type=truvari.restricted_int, default=50,
                        help="Minimum size of entry to annotate (%(default)s)")
    truvari.setup_logging()
    return parser.parse_args(args)


def edit_header(my_vcf):
    """
    Add INFO for new fields to vcf
    """
    header = my_vcf.header.copy()
    header.add_line(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SVTYPE">')
    header.add_line(
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SVLEN">')
    return header


def svinfo_main(cmdargs):
    """
    Main method
    """
    args = parse_args(cmdargs)
    vcf = pysam.VariantFile(args.input)
    n_header = edit_header(vcf)
    with pysam.VariantFile(args.output, 'w', header=n_header) as out:
        for entry in vcf:
            sz = truvari.entry_size(entry)
            if sz >= args.minsize:
                entry.translate(n_header)
                svtype = truvari.entry_variant_type(entry)
                entry.info["SVTYPE"] = svtype
                entry.info["SVLEN"] = sz
            out.write(entry)
    logging.info("Finished svinfo")
