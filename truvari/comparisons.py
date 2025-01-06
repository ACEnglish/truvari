"""
Collection of methods with event helpers
that compare events, coordinates, or transform vcf entries
"""
import edlib


def coords_within(qstart, qend, rstart, rend, end_within):
    """
    Returns if a span is within the provided [start, end). All coordinates assumed 0-based
    One exception is made for REPL types with anchor bases. This is for TR callers and GIABTR benchmark
    where variants will span the entire includebed region.

    :param `qstart`: query start position
    :type `qstart`: integer
    :param `qend`: query end position
    :type `qend`: integer
    :param `start`: start of span
    :type `start`: integer
    :param `end`: end of span
    :type `end`: integer
    :param `end_within`:  if true, qend <= rend, else qend < rend
    :type `end_within`: bool

    :return: If the coordinates are within the span
    :rtype: bool
    """
    # REPL types, probably
    if qstart in (rstart - 1, rstart) and qend == rend:
        return True
    if end_within:
        ending = qend <= rend
    else:
        ending = qend < rend
    return qstart >= rstart and ending


def best_seqsim(a_seq, b_seq, st_dist):
    """
    Returns best of roll, unroll, and direct sequence similarity
    """
    return max(roll_seqsim(a_seq, b_seq), unroll_seqsim(a_seq, b_seq, st_dist),
               unroll_seqsim(b_seq, a_seq, -st_dist), seqsim(a_seq, b_seq))


def roll_seqsim(a_seq, b_seq):
    """
    Compare the lexicographically smallest rotations of two sequences
    """
    def smallest_rotation(s):
        doubled = s + s
        n = len(s)
        return min(doubled[i:i + n] for i in range(n))
    r_a = smallest_rotation(a_seq)
    r_b = smallest_rotation(b_seq)
    return seqsim(r_a, r_b)


def overlap_percent(astart, aend, bstart, bend):
    """
    Calculates the percent of range A which overlaps with range B

    :param `astart`: First range's start position
    :type `astart`: int
    :param `aend`: First range's end position
    :type `aend`: int
    :param `bstart`: Second range's start position
    :type `bstart`: int
    :param `bend`: Second range's end position
    :type `bend`: int

    :return: overlap percent
    :rtype: float
    """
    if astart >= bstart and aend <= bend:
        return 1
    ovl_start = max(astart, bstart)
    ovl_end = min(aend, bend)
    if ovl_start < ovl_end:  # Otherwise, they're not overlapping
        ovl_pct = float(ovl_end - ovl_start) / (aend - astart)
    else:
        ovl_pct = 0
    return ovl_pct


def overlaps(s1, e1, s2, e2):
    """
    Check if two start/end ranges have overlap

    :param `s1`: range 1 start
    :type `s1`: int
    :param `e1`: range 1 end
    :type `e1`: int
    :param `s2`: range 2 start
    :type `s2`: int
    :param `e2`: range 2 end
    :type `e2`: int

    :return: True if ranges overlap
    :rtype: bool
    """
    s_cand = max(s1, s2)
    e_cand = min(e1, e2)
    return s_cand < e_cand


def reciprocal_overlap(astart, aend, bstart, bend):
    """
    Calculates reciprocal overlap of two ranges

    :param `astart`: First range's start position
    :type `astart`: int
    :param `aend`: First range's end position
    :type `aend`: int
    :param `bstart`: Second range's start position
    :type `bstart`: int
    :param `bend`: Second range's end position
    :type `bend`: int

    :return: reciprocal overlap
    :rtype: float
    """
    ovl_start = max(astart, bstart)
    ovl_end = min(aend, bend)
    if ovl_start < ovl_end:  # Otherwise, they're not overlapping
        ovl_pct = float(ovl_end - ovl_start) / \
            max(aend - astart, bend - bstart)
    else:
        ovl_pct = 0
    return ovl_pct


def seqsim(allele1, allele2):
    """
    Calculate similarity of two sequences

    :param `allele1`: first entry
    :type `allele1`: str
    :param `allele2`: second entry
    :type `allele2`: str

    :return: sequence similarity
    :rtype: float
    """
    allele1 = allele1.upper()
    allele2 = allele2.upper()
    scr = edlib.align(allele1, allele2)
    totlen = len(allele1) + len(allele2)
    return (totlen - scr["editDistance"]) / totlen


def sizesim(sizeA, sizeB):
    """
    Calculate the size similarity percent and size diff for two sizes

    :param `sizeA`: first size
    :type `sizeA`: int
    :param `sizeB`: second size
    :type `sizeB`: int

    :return: size similarity percent and size diff (A - B)
    :rtype: (float, int)
    """
    if sizeA == 0 or sizeB == 0:
        if sizeA == sizeB:
            return 1, 0
        sizeA = max(sizeA, 1)
        sizeB = max(sizeB, 1)
    return min(sizeA, sizeB) / float(max(sizeA, sizeB)), sizeA - sizeB


def unroll_seqsim(seqA, seqB, p):
    """
    Unroll two sequences and compare.
    See https://gist.github.com/ACEnglish/1e7421c46ee10c71bee4c03982e5df6c for details

    :param `seqA`: sequence upstream sequence
    :type `seqA`: string
    :param `seqB`: sequence downstream sequence
    :type `seqB`: string
    :param `p`: how many bases upstream seqA is from seqB
    :type `p`: integer

    :return: sequence similarity of seqA vs seqB after unrolling
    :rtype: float
    """
    f = p % len(seqB)
    uB = seqB[-f:] + seqB[:-f]
    return seqsim(seqA, uB)
