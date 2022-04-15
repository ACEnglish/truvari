"""
Wrapper around the different annotations available
"""
import argparse
import truvari.annos as tannos

ANNOS = {"gcpct": tannos.gcpct_main,
         "gtcnt": tannos.gtcnt_main,
         "trf": tannos.trf_main,
         "grm": tannos.grm_main,
         "repmask": tannos.rmk_main,
         "remap": tannos.remap_main,
         "hompct": tannos.hompct_main,
         "numneigh": tannos.numneigh_main,
         "svinfo": tannos.svinfo_main,
         "bpovl": tannos.bpovl_main,
         "density": tannos.density_main,
         "dpcnt": tannos.dpcnt_main,
         "lcr": tannos.lcr_main}

USAGE = """\
Truvari annotations:
        gcpct, gtcnt, trf, grm, repmask, remap, hompct, numneigh, svinfo, bpovl, density, dpcnt, lcr
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
