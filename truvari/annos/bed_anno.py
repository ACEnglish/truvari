"""
Annotate over a generic bedfile

Need to build the 'tree maker' first?

Need to make sure and set it so
    --multimatch

And also, need to specify that col[3] (bed name)
must be the INFO=;oaiwef
and the header lines "^#" must be the header information
"""
import argparse
import truvari
from acebinf import setup_logging


def main():
    """
    Main
    """
    self.refanno = truvari.make_bedanno_tree(refanno) if refanno else None

    srep_hits = self.refanno[0][entry.chrom].overlap(entry.start, entry.stop)
    for i in srep_hits:
        i = i.data
        lookup[i["SREP_repeats"][0]] = i["SREP_copies"][0]

    if srep_hits is not None:
        for i in srep_hits:
            for k,v in i.data.items():
                n_dat[k].extend(v)

def parse_args(args):
    """
    Pull the command line parameters
    """
    def restricted_float(x):
        x = float(x)
        if x < 0.0 or x > 1.0:
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
        return x
    parser = argparse.ArgumentParser(prog="trf", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, default="/dev/stdin",
                        help="VCF to annotate (stdin)")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output filename (stdout)")
    parser.add_argument("-e", "--executable", type=str, default="trf409.linux64",
                        help="Path to tandem repeat finder (%(default)s)")
    parser.add_argument("-f", "--full", action="store_true",
                        help="Write full trf output to entries")
    parser.add_argument("-m", "--min-length", type=int, default=50,
                        help="Minimum size of entry to annotate (%(default)s)")
    parser.add_argument("-t", "--threshold", type=restricted_float, default=.8,
                        help="Threshold for pct of allele covered (%(default)s)")
    parser.add_argument("-T", "--trf-params", type=str, default="2 7 7 80 10 50 500 -m -f -h -d -ngs",
                        help="Default parameters to send to trf (%(default)s)")
    parser.add_argument("-R", "--ref-bed", type=str, default=None,
                        help="Reference bed of simple repeats")
    parser.add_argument("-r", "--ref", type=str, default=None,
                        help="Reference fasta file (use with -R)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    setup_logging(args.debug)
    return args


