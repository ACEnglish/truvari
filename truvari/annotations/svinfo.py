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
    parser.add_argument("input", nargs="?", type=str, default="/dev/stdin",
                        help="VCF to annotate (stdin)")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output filename (stdout)")
    parser.add_argument("-m", "--minsize", type=truvari.restricted_int, default=50,
                        help="Minimum size of entry to annotate (%(default)s)")
    truvari.setup_logging(show_version=True)
    return parser.parse_args(args)


def edit_header(header):
    """
    Add INFO for new fields to vcf
    """
    header.add_line(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SVTYPE">')
    header.add_line(
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SVLEN">')
    return header


def add_svinfo(entry, min_size=0, n_header=None):
    """
    Add svinfo
    """
    if "SVTYPE" in entry.info:
        del entry.info['SVTYPE']
    if "SVLEN" in entry.info:
        del entry.info['SVLEN']
    sz = truvari.entry_size(entry)
    if sz < min_size:
        return
    if n_header:
        entry.translate(n_header)
    svtype = truvari.entry_variant_type(entry)
    entry.info["SVTYPE"] = svtype.name
    entry.info["SVLEN"] = sz


def svinfo_main(cmdargs):
    """
    Main method
    """
    args = parse_args(cmdargs)
    vcf = pysam.VariantFile(args.input)
    n_header = edit_header(vcf.header.copy())
    with pysam.VariantFile(args.output, 'w', header=n_header) as out:
        for entry in vcf:
            add_svinfo(entry, args.minsize, n_header)
            out.write(entry)
    logging.info("Finished svinfo")
