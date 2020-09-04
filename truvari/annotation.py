"""
Wrapper around the different annotations available
"""
import argparse
from truvari.annos import *

ANNOS = {"gcpct": gcpct_main,
         "gtcnt": gtcnt_main,
         "trf": trf_main,
         "grm": grm_main,
         "repmask": rmk_main,
         "remap": remap_main,
         "hompct": hompct_main}

USAGE = """\
Truvari annotations:
        gcpct, gtcnt, trf, grm, repmask, remap, hompct
"""

def parseArgs(args):
    """
    Argument parsing
    """
    parser = argparse.ArgumentParser(prog="truvari anno", description=USAGE,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("cmd", metavar="CMD", choices=ANNOS.keys(), type=str,
                        help="Annotation to run")
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help="Options to pass to the annotation")


    args = parser.parse_args(args)
    return args

def anno_main(args):
    """
    Simple wrapper around the main
    """
    args = parseArgs(args)
    ANNOS[args.cmd](args.options)
