"""
VCF Stats for SVs with sizebins and svtypes
"""
import sys
import argparse
import warnings

from collections import defaultdict, Counter

import vpq
import numpy
import pysam
import joblib
import truvari


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
    # need to add qualbin scaling
    return parser.parse_args(args)


def get_svtype(svtype):
    """
    Turn to a vpq svtype
    """
    try:
        return eval(f"vpq.SV.{svtype}")
    except AttributeError:
        pass
    return vpq.SV.UNK


def get_sizebin(sz):
    """
    Bin a given size
    """
    sz = abs(sz)
    for key, maxval in zip(vpq.SZBINS, vpq.SZBINMAX):
        if sz <= maxval:
            return key
    return None


def get_gt(gt):
    """
    return GT enum
    """
    if None in gt:
        return vpq.GT.NON
    if len(gt) != 2:
        return vpq.GT.UNK
    if gt == (0, 0):
        return vpq.GT.REF
    if gt[0] != gt[1]:
        return vpq.GT.HET
    if gt[0] == gt[1]:
        return vpq.GT.HOM
    return vpq.GT.UNK


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


def stats_main(cmdargs):
    """
    Main stats
    """
    warnings.filterwarnings('ignore')
    args = parse_args(cmdargs)

    output = open(args.out, 'w')

    data = numpy.zeros((len(vpq.SV), len(vpq.SZBINS), len(vpq.GT), len(vpq.QUALBINS)))

    for fn in args.VCF:
        vcf = pysam.VariantFile(fn)
        for entry in vcf:
            sv = get_svtype(truvari.entry_variant_type(entry))
            sz = get_sizebin(truvari.entry_size(entry))
            gt = get_gt(entry.samples[0]["GT"])
            qual, idx = get_scalebin(entry.qual)
            data[sv.value, vpq.SZBINS.index(sz), gt.value, idx] += 1

    if args.dataframe:
        joblib.dump(data, args.dataframe)
    output.write("# SV counts\n")
    for i in vpq.SV:
        output.write("%s\t%d\n" % (i.name, data[i.value].sum()))

    output.write("# SVxSZ counts\n")
    output.write("\t%s\n" % "\t".join(x.name for x in vpq.SV))
    for idx, sz in enumerate(vpq.SZBINS):
        output.write(sz)
        for sv in vpq.SV:
            output.write("\t%d\n" % (data[sv.value, idx].sum()))
        output.write('\n')

    output.write("# GT counts\n")
    for i in vpq.GT:
        output.write("%s\t%d\n" % (i.name, data[:, :, i.value].sum()))

    output.write("# SVxGT counts\n")
    output.write("\t%s\n" % "\t".join(x.name for x in vpq.SV))
    for gt in vpq.GT:
        output.write(gt.name)
        for sv in vpq.SV:
            output.write("\t%d" % (data[sv.value, :, gt.value].sum()))
        output.write('\n')

    output.write("# Het/Hom ratio\n")
    output.write("%f\n" % (data[:, :, vpq.GT.HET.value].sum() / data[:, :, vpq.GT.HOM.value].sum()))

    output.write("# SVxSZ Het/Hom ratios\n")
    output.write("\t%s\n" % "\t".join(x.name for x in vpq.SV))
    for idx, sz in enumerate(vpq.SZBINS):
        output.write(sz)
        for sv in vpq.SV:
            het = data[sv.value, idx, vpq.GT.HET.value].sum()
            hom = data[sv.value, idx, vpq.GT.HOM.value].sum()
            output.write("\t%.2f" % (het / hom))
        output.write("\n")

    output.write("# QUAL distribution\n")
    for pos, i in enumerate(range(0, 100, 10)):
        output.write("[%d,%d)\t%d\n" % (i, i+10, data[:, :, :, pos].sum()))
    output.write(">=%d\t%d\n" % (100, data[:, :, :, pos + 1].sum()))


if __name__ == '__main__':
    main()
