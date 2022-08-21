"""
Wrapper around the different annotations available
"""
import argparse
from rich.console import Console
import truvari.annotations as tannos

ANNOS = {"gcpct": ("GC Percent", tannos.gcpct_main),
         "gtcnt": ("Genotype Counts", tannos.gtcnt_main),
         "trf": ("Tandem Repeats", tannos.trf_main),
         "grm": ("Mappability", tannos.grm_main),
         "repmask": ("Repeats", tannos.rmk_main),
         "remap": ("Allele Remapping", tannos.remap_main),
         "hompct": ("Homozygous Percent", tannos.hompct_main),
         "numneigh": ("Number of Neighbors", tannos.numneigh_main),
         "svinfo": ("SVINFO Fields", tannos.svinfo_main),
         "bpovl": ("Annotation Intersection", tannos.bpovl_main),
         "density": ("Variant Density", tannos.density_main),
         "dpcnt": ("Call Depth Counts", tannos.dpcnt_main),
         "lcr": ("Low-complexity Regions", tannos.lcr_main),
         "grpaf": ("Sample Group Allele Frequency", tannos.grpaf_main)}


USAGE = "Truvari annotations:\n" + "\n".join([f"    [bold][cyan]{k:9}[/][/] {t[0]}" for k,t in ANNOS.items()])

class ArgumentParser(argparse.ArgumentParser):
    """
    Custom ArgumentParser
    """
    def _print_message(self, message, file=None):
        """ pretty print """
        console = Console(stderr=True)
        console.print(message, highlight=False)

def parseArgs(args):
    """
    Argument parsing
    """
    parser = ArgumentParser(prog="truvari anno", description=USAGE,
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
    ANNOS[args.cmd][1](args.options)
