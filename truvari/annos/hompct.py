"""
Calcluate the the Hom / (Het + Hom) of variants in the region of SVs
Requires the VCF to contain SVs beside SNPs/Indels
"""
import logging
import argparse

import pysam
import truvari


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
    parser.add_argument("-b", "--buffer", type=truvari.restricted_int, default=5000,
                        help="Number of base-pairs up/dn-stream to query (%(default)s)")
    parser.add_argument("-m", "--minanno", type=truvari.restricted_int, default=50,
                        help="Minimum size of event to annotate (%(default)s)")
    parser.add_argument("-M", "--maxgt", type=truvari.restricted_int, default=1,
                        help="Largest event size to count for genotyping (%(default)s)")
    parser.add_argument("-c", "--mincount", type=truvari.restricted_int, default=0,
                        help="Minimum number of genotyping events to report HOMPCT (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    return args


def hompct_main(cmd_args):
    """
    Main
    """
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
        if tot < args.mincount:
            return None

        if tot == 0:
            return float('nan')

        return float(format((homs / tot) * 100, ".1f"))

    header = v.header.copy()
    header.add_line(('##INFO=<ID=HOMPCT,Number=1,Type=Float,'  # pylint: disable=consider-using-f-string
                     'Description="Percent of calls < %dbp long within %dbp that are homozygous">')
                    % (args.maxgt, args.buffer))

    out = pysam.VariantFile(args.output, 'w', header=header)
    v2 = pysam.VariantFile(args.input)
    for entry in v2:
        if truvari.entry_size(entry) >= args.minanno:
            entry.translate(header)
            anno = get_pct(entry.chrom, *truvari.entry_boundaries(entry))
            if anno is not None:
                entry.info["HOMPCT"] = anno
        out.write(entry)
    logging.info("Finished hompct")
