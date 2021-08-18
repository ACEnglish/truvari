#!/usr/bin/env python
"""
Truvari main entrypoint
"""
import sys
import argparse

from truvari import __version__
from truvari.bench import bench_main
from truvari.collapse import collapse_main
from truvari.consistency_report import consistency_main
from truvari.annotation import anno_main
from truvari.vcf2df import vcf2df_main

# pylint: disable=unused-argument
def in_progress(args):
    """placeholder"""
    print('working on it...')

def version(args):
    """Print the version"""
    print("Truvari v%s" % __version__)

TOOLS = {'bench': bench_main,
         'consistency': consistency_main,
         'anno': anno_main,
         'collapse': collapse_main,
         'vcf2df': vcf2df_main,
         'version': version}

# create-anno      Create an index of an annotation source for reuse

USAGE = """\
Truvari v%s - Structural Variant Benchmarking and Annotation

    CMDs:
        bench         Performance metrics from comparison of two VCFs
        consistency   Consistency report between multiple VCFs
        anno          Annotate a VCF
        collapse      Collapse possibly redundant VCF entries
        vcf2df        Turn a VCF into a pandas DataFrame
        version       Print the Truvari version and exit
""" % __version__


def main():
    """
    Argument parsing
    """
    parser = argparse.ArgumentParser(prog="truvari", description=USAGE,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("cmd", metavar="CMD", choices=TOOLS.keys(), type=str,
                        help="Command to execute")
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help="Options to pass to the command")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()

    args = parser.parse_args()

    TOOLS[args.cmd](args.options)

if __name__ == '__main__':
    main()
