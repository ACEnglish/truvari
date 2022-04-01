#!/usr/bin/env python3
"""
Truvari main entrypoint
"""
import sys
import argparse

import truvari
from truvari import __version__
from truvari.bench import bench_main
from truvari.divide import divide_main
from truvari.vcf2df import vcf2df_main
from truvari.annotation import anno_main
from truvari.collapse import collapse_main
from truvari.segmentation import segment_main
from truvari.consistency_report import consistency_main

# pylint: disable=unused-argument


def version(args):
    """Print the version"""
    print(f"Truvari v{__version__}")


TOOLS = {'bench': bench_main,
         'consistency': consistency_main,
         'anno': anno_main,
         'collapse': collapse_main,
         'vcf2df': vcf2df_main,
         'segment': segment_main,
         'divide': divide_main,
         'version': version}

USAGE = f"""\
Truvari v{__version__} - Structural Variant Benchmarking and Annotation

    CMDs:
        bench         Performance metrics from comparison of two VCFs
        consistency   Consistency report between multiple VCFs
        anno          Annotate a VCF
        collapse      Collapse possibly redundant VCF entries
        vcf2df        Turn a VCF into a pandas DataFrame
        segment       Normalization of SVs into disjointed genomic regions
        divide        Divide a VCF into multiple parts
        version       Print the Truvari version and exit
"""


class ArgumentParser(argparse.ArgumentParser):
    """
    Custom argument parser error
    """

    def error(self, message):
        """
        Check for similar commands before exiting
        """
        sys.stderr.write(f'{self.prog}: error: {message}\n')
        if message.startswith("argument CMD: invalid choice"):
            guess = truvari.help_unknown_cmd(sys.argv[1], TOOLS.keys())
            if guess:
                sys.stderr.write(f"\nThe most similar command is\n\t{guess}\n")

        self.exit(2)


def main():
    """
    Main entrypoint for truvari tools
    """
    parser = ArgumentParser(prog="truvari", description=USAGE,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("cmd", metavar="CMD", choices=TOOLS.keys(), type=str, default=None,
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
