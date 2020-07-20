"""
VCF Stats for SVs with sizebins and svtypes
"""
import sys
import argparse
import warnings

from enum import Enum
from collections import defaultdict, Counter

import numpy
import pysam
import joblib
import truvari

SZBINS = ["[0,50)", "[50,100)", "[100,200)", "[200,300)", "[300,400)",
          "[400,600)", "[600,800)", "[800,1k)", "[1k,2.5k)",
          "[2.5k,5k)", ">=5k"]
SZBINMAX = [50, 100, 200, 300, 400, 600, 800, 1000, 2500, 5000, sys.maxsize]
QUALBINS = [f"[{x},{x+10})" for x in range(0, 100, 10)] + [">=100"]

class GT(Enum):
    """ Genotypes """
    NON = 3
    REF = 0
    HET = 1
    HOM = 2
    UNK = 4


class SV(Enum):
    """ SVtypes """
    DEL = 0
    INS = 1
    DUP = 2
    INV = 3
    NON = 4  # Not and SV, SVTYPE
    UNK = 5  # Unknown SVTYPE


def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="stats", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("VCF", nargs="+", default="/dev/stdin",
                        help="VCFs to annotate (stdin)")
    parser.add_argument("-d", "--dataframe", default=None,
                        help="Write dataframe joblib to file (%(default)s)")
    parser.add_argument("-o", "--out", default="/dev/stdout",
                        help="Stats output file (%(default)s")
    parser.add_argument("--qmin", default=0, type=int,
                        help="Minimum QUAL score found in VCF (%(default)s)")
    parser.add_argument("--qmax", default=None, type=int,
                        help="Maximum QUAL score found in VCF (None)")

    args = parser.parse_args(args)
    if args.qmax is None:
        args.qmax = sys.maxsize
    return args


def get_svtype(svtype):
    """
    Turn to a svtype
    """
    try:
        return eval(f"SV.{svtype}")
    except AttributeError:
        pass
    return SV.UNK


def get_sizebin(sz):
    """
    Bin a given size
    """
    sz = abs(sz)
    for key, maxval in zip(SZBINS, SZBINMAX):
        if sz <= maxval:
            return key
    return None


def get_gt(gt):
    """
    return GT enum
    """
    if None in gt:
        return GT.NON
    if len(gt) != 2:
        return GT.UNK
    if gt == (0, 0):
        return GT.REF
    if gt[0] != gt[1]:
        return GT.HET
    if gt[0] == gt[1]:
        return GT.HOM
    return GT.UNK


def get_scalebin(x, rmin=0, rmax=100, tmin=0, tmax=100, step=10):
    """
    Scale variable x from rdomain to tdomain with step sizes
    return key, index

    rmin denote the minimum of the range of your measurement
    tmax denote the maximum of the range of your measurement
    tmin denote the minimum of the range of your desired target scaling
    tmax denote the maximum of the range of your desired target scaling
    """
    newx = (x - rmin) / (rmax - rmin) * (tmax - tmin) + tmin
    pos = 0
    for pos, i in enumerate(range(tmin, tmax, step)):
        if newx < i + step:
            return "[%d,%d)" % (i, i+step), pos
    return ">=%d" % (tmax), pos + 1

def format_stats(data, sample=False, output=sys.stdout):
    """
    Given a stats array, output
    """
    output.write("# SV counts\n")
    for i in SV:
        if sample:
            output.write("%s\t%d\n" % (i.name, data[i.value, :, :, GT.HET.value:GT.HOM.value + 1].sum()))
        else:
            output.write("%s\t%d\n" % (i.name, data[i.value].sum()))
    if sample:
        output.write("%s\t%d\n" % ("total", data[:, :, :, GT.HET.value:GT.HOM.value + 1].sum()))
    else:
        output.write("%s\t%d\n" % ("total", data.sum()))


    output.write("# SVxSZ counts\n")
    output.write("\t%s\n" % "\t".join(x.name for x in SV))
    for idx, sz in enumerate(SZBINS):
        output.write(sz)
        for sv in SV:
            if sample:
                output.write("\t%d" % (data[sv.value, idx, :, GT.HET.value:GT.HOM.value + 1].sum()))
            else:
                output.write("\t%d" % (data[sv.value, idx].sum()))
        output.write('\n')

    output.write("# QUAL distribution\n")
    for pos, i in enumerate(range(0, 100, 10)):
        if sample:
            output.write("[%d,%d)\t%d\n" % (i, i+10, data[:, :, pos, GT.HET.value:GT.HOM.value + 1].sum()))
        else:
            output.write("[%d,%d)\t%d\n" % (i, i+10, data[:, :, pos].sum()))
    output.write(">=%d\t%d\n" % (100, data[:, :, pos + 1].sum()))


    if sample:
        output.write("# GT counts\n")
        for i in GT:
            output.write("%s\t%d\n" % (i.name, data[:, :, :, i.value].sum()))

        output.write("# SVxGT counts\n")
        output.write("\t%s\n" % "\t".join(x.name for x in SV))
        for gt in GT:
            output.write(gt.name)
            for sv in SV:
                output.write("\t%d" % (data[sv.value, :, :, gt.value].sum()))
            output.write('\n')

        output.write("# Het/Hom ratio\n")
        output.write("%f\n" % (data[:, :, :, GT.HET.value].sum() / data[:, :, :, GT.HOM.value].sum()))

        output.write("# SVxSZ Het/Hom ratios\n")
        output.write("\t%s\n" % "\t".join(x.name for x in SV))

        for idx, sz in enumerate(SZBINS):
            output.write(sz)
            for sv in SV:
                het = data[sv.value, idx, :, GT.HET.value].sum()
                hom = data[sv.value, idx, :, GT.HOM.value].sum()
                output.write("\t%.2f" % (het / hom))
            output.write("\n")

def generate_stat_table(vcf_fn, args):
    """
    Given a vcf filename, create a numpy array with dimensions counting
    [SVTYPE, SZBINS, GT, QUALBINS]
    """

    vcf = pysam.VariantFile(vcf_fn)
    ret = {}
    for i in vcf.header.samples:
        ret[i] = numpy.zeros((len(SV), len(SZBINS), len(QUALBINS), len(GT)))
    ret["total"] = numpy.zeros((len(SV), len(SZBINS), len(QUALBINS)))
    for entry in vcf:
        sv = get_svtype(truvari.entry_variant_type(entry))
        sz = get_sizebin(truvari.entry_size(entry))
        if entry.qual is not None:
            qual, idx = get_scalebin(entry.qual, args.qmin, args.qmax)
        else:
            qual, idx = 0, 0
        for i in vcf.header.samples:
            gt = get_gt(entry.samples[i]["GT"])
            ret[i][sv.value, SZBINS.index(sz), idx, gt.value] += 1
        ret["total"][sv.value, SZBINS.index(sz), idx] += 1

    return ret


def stats_main(cmdargs):
    """
    Main stats
    """
    warnings.filterwarnings('ignore')
    args = parse_args(cmdargs)

    output = open(args.out, 'w')
    data = None
    for vcf in args.VCF:
        cur = generate_stat_table(vcf, args)
        if data is None:
            data = cur
        else:
            for i in cur:
                data[i] += cur[i]

    output.write("## Total Stats:\n")
    format_stats(data["total"], False, output)

    for i in [_ for _ in data if _ != 'total']:
        output.write("## %s Stats:\n" % (i))
        format_stats(data[i], True, output)

    if args.dataframe:
        joblib.dump(data, args.dataframe)

if __name__ == '__main__':
    main()
