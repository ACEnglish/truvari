"""
VCF Stats for SVs with sizebins and svtypes
"""
import sys

from enum import Enum

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
    NON = 4  # Not an SV, SVTYPE
    UNK = 5  # Unknown SVTYPE

def get_svtype(svtype):
    """
    Turn to a svtype
    """
    try:
        return truvari.SV.__members__[svtype]
    except AttributeError:
        pass
    return SV.UNK


def get_sizebin(sz):
    """
    Bin a given size
    """
    sz = abs(sz)
    for key, maxval in zip(SZBINS, SZBINMAX):
        if sz < maxval:
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
