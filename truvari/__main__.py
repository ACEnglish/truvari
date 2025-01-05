#!/usr/bin/env python3
"""
Truvari main entrypoint
"""
import sys
import argparse
from importlib.metadata import version

from rich.console import Console

import truvari
from truvari import __version__
from truvari.anno import anno_main
from truvari.phab import phab_main
from truvari.bench import bench_main
from truvari.divide import divide_main
from truvari.vcf2df import vcf2df_main
from truvari.refine import refine_main
from truvari.collapse import collapse_main
from truvari.stratify import stratify_main
from truvari.segmentation import segment_main
from truvari.consistency import consistency_main
from truvari.make_ga4gh import make_ga4gh_main

def flat_version(args):
    """Print the version"""
    if len(args) and args[0].count("-v"):
        print(f"Truvari {version('truvari')}")
    else:
        print(f"Truvari v{__version__}")


TOOLS = {"bench": bench_main,
         "consistency": consistency_main,
         "anno": anno_main,
         "collapse": collapse_main,
         "vcf2df": vcf2df_main,
         "segment": segment_main,
         "stratify": stratify_main,
         "divide": divide_main,
         "phab": phab_main,
         "refine": refine_main,
         "ga4gh": make_ga4gh_main,
         "version": flat_version}

USAGE = f"""\
[bold]Truvari v{__version__}[/] Structural Variant Benchmarking and Annotation

Available commands:
    [bold][cyan]bench[/][/]         Performance metrics from comparison of two VCFs
    [bold][cyan]collapse[/][/]      Collapse possibly redundant VCF entries
    [bold][cyan]anno[/][/]          Annotate a VCF
    [bold][cyan]consistency[/][/]   Consistency report between multiple VCFs
    [bold][cyan]vcf2df[/][/]        Turn a VCF into a pandas DataFrame
    [bold][cyan]segment[/][/]       Normalization of SVs into disjointed genomic regions
    [bold][cyan]stratify[/][/]      Count variants per-region in vcf
    [bold][cyan]divide[/][/]        Divide a VCF into independent shards
    [bold][cyan]phab[/][/]          Variant harmonization using MSA
    [bold][cyan]refine[/][/]        Automated bench result refinement with phab
    [bold][cyan]ga4gh[/][/]         Convert Truvari result to GA4GH
    [bold][cyan]version[/][/]       Print the Truvari version and exit
"""


class ArgumentParser(argparse.ArgumentParser):
    """
    Custom argument parser error
    """

    def error(self, message):
        """
        Check for similar commands before exiting
        """
        console = Console(stderr=True)
        console.print(f'{self.prog}: error: {message}')
        if message.startswith("argument CMD: invalid choice"):
            guess = truvari.help_unknown_cmd(sys.argv[1], TOOLS.keys())
            if guess:
                console.print(f"\nThe most similar command is\n\t[bold][cyan]{guess}[/][/]\n")

        self.exit(2)

    def _print_message(self, message, file=None):
        """ pretty print """
        console = Console(stderr=True)
        console.print(message, highlight=False)

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
