"""
Collection of methods with event helpers
that compare events, coordinates, or transform vcf entries
"""
import re
import hashlib

import edlib

import truvari


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


def entry_is_present(entry, sample=None, allow_missing=True):
    """
    Checks if entry's sample genotype is present and is heterozygous or homozygous (a.k.a. present)
    If allow_missing, just check for a 1 in the genotype. Otherwise, a missing ('.') genotype isn't
    considered present

    :param `entry`: entry to check
    :type `entry`: :class:`pysam.VariantRecord`
    :param `sample`: sample name
    :type `sample`: string, optional

    :return: True if variant is present in the sample
    :rtype: bool

    Example
        >>> import truvari
        >>> import pysam
        >>> v = pysam.VariantFile('repo_utils/test_files/variants/input1.vcf.gz')
        >>> truvari.entry_is_present(next(v))
        True
    """
    if sample is None:
        sample = entry.samples.keys()[0]
    if allow_missing:
        return 1 in entry.samples[sample]["GT"]
    return "GT" in entry.samples[sample] and \
           truvari.get_gt(entry.samples[sample]["GT"]) in [
        truvari.GT.HET, truvari.GT.HOM]


def entry_seq_similarity(entryA, entryB, roll=True):
    """
    Calculate sequence similarity of two entries. If reference is not None,
    compare their shared reference context. Otherwise, use the unroll technique.

    :param `entryA`: first entry
    :type `entryA`: :class:`pysam.VariantRecord`
    :param `entryB`: second entry
    :type `entryB`: :class:`pysam.VariantRecord`
    :param `roll`: compare the lexicographically minimum rotation of sequences
    :type `roll`: bool

    :return: sequence similarity
    :rtype: float
    """
    # Shortcut to save compute - probably unneeded
    if entryA.ref == entryB.ref and entryA.alts[0] == entryB.alts[0]:
        return 1.0

    # Inversions aren't rolled
    if (entryA.svtype() == truvari.SV.INV and entryB.svtype() == truvari.SV.INV):
        allele1 = entryA.alts[0]
        allele2 = entryB.alts[0]
        return seqsim(allele1, allele2)

    a_seq = entryA.ref if entryA.svtype() == truvari.SV.DEL else entryA.alts[0]
    a_seq = a_seq.upper()
    b_seq = entryB.ref if entryB.svtype() == truvari.SV.DEL else entryB.alts[0]
    b_seq = b_seq.upper()
    st_dist, ed_dist = entryA.distance(entryB)

    if not roll or st_dist == 0 or ed_dist == 0:
        return seqsim(a_seq, b_seq)

    if st_dist < 0:
        st_dist *= -1
    else:
        a_seq, b_seq = b_seq, a_seq

    # Return best of rolled, unrolled from both ends, and direct similarity
    # Whichever is highest is how similar these sequences can be
    return best_seqsim(a_seq, b_seq, st_dist)


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


def entry_reciprocal_overlap(entry1, entry2, ins_inflate=True):
    """
    Calculates reciprocal overlap of two entries

    :param `entry1`: First entry
    :type `entry1`: :class:`pysam.VariantRecord`
    :param `entry2`: Second entry
    :type `entry2`: :class:`pysam.VariantRecord`
    :param `ins_inflate`: inflate entry boundaries
    :type `ins_inflate`: bool

    :return: The reciprocal overlap
    :rtype: float

    Example
        >>> import truvari
        >>> import pysam
        >>> v = pysam.VariantFile('repo_utils/test_files/variants/input1.vcf.gz')
        >>> a = next(v)
        >>> b = next(v)
        >>> truvari.entry_reciprocal_overlap(a, b)
        0
    """
    astart, aend = entry1.boundaries(ins_inflate)
    bstart, bend = entry2.boundaries(ins_inflate)
    return reciprocal_overlap(astart, aend, bstart, bend)


def entry_same_variant_type(entryA, entryB, dup_to_ins=False):
    """
    Check if entryA svtype == entryB svtype

    :param `entryA`: first entry
    :type `entryA`: :class:`pysam.VariantRecord`
    :param `entryB`: second entry
    :type `entryB`: :class:`pysam.VariantRecord`
    :param `dup_to_ins`: Convert DUP to INS types
    :type `dup_to_ins`: bool

    :return: True if entry SVTYPEs match
    :rtype: bool
    """
    a_type = entryA.svtype()
    b_type = entryB.svtype()
    if dup_to_ins and a_type == truvari.SV.DUP:
        a_type = truvari.SV.INS
    if dup_to_ins and b_type == truvari.SV.DUP:
        b_type = truvari.SV.INS
    return a_type == b_type


def entry_size_similarity(entryA, entryB):
    """
    Calculate the size similarity and difference for two entries

    :param `entryA`: first entry
    :type `entryA`: :class:`pysam.VariantRecord`
    :param `entryB`: second entry
    :type `entryB`: :class:`pysam.VariantRecord`

    :return: size similarity and size diff (A - B)
    :rtype: (float, int)

    Example
        >>> import truvari
        >>> import pysam
        >>> v = pysam.VariantFile('repo_utils/test_files/variants/input1.vcf.gz')
        >>> a = next(v)
        >>> b = next(v)
        >>> truvari.entry_size_similarity(a, b)
        (0.07142857142857142, 13)
    """
    sizeA = entry_size(entryA)
    sizeB = entry_size(entryB)
    return sizesim(sizeA, sizeB)


def entry_to_hash(entry, hasher=hashlib.sha1):
    """
    Turn variant into a key and hash with provided hasher

    :param `entry`: entry
    :type `entry`: :class:`pysam.VariantRecord`
    :param `hasher`: hashing function
    :type `entry`: method

    :return: hash
    :rtype: string
    """
    return hasher(entry_to_key(entry).encode()).hexdigest()


def entry_to_key(entry, prefix="", bounds=False):
    """
    Turn a vcf entry into a hashable key string. Use the prefix (base/comp) to indicate the source
    VCF when consolidating multiple files' calls. If bounds: call entry_boundaries for start/stop.

    .. warning::
        If a caller redundantly calls a variant exactly the same it will not have a unique key

    :param `entry`: entry to stringify
    :type `entry`: :class:`pysam.VariantRecord`
    :param `prefix`: prefix
    :type `prefix`: string, optional
    :param `bounds`: use entry_boundaries
    :type `bounds`: bool, optional

    :return: hashable string uniquely identifying the variant
    :rtype: string
    """
    if prefix:
        prefix += '.'
    alt = entry.alts[0] if entry.alts is not None else "."
    if bounds:
        start, end = entry.boundaries()
        return f"{prefix}{entry.chrom}:{start}-{end}({entry.ref}|{alt})"
    return f"{prefix}{entry.chrom}:{entry.start}-{entry.stop}.{alt}"


def entry_within_tree(entry, tree):
    """
    Extract entry and tree boundaries to call `coords_within`
    """
    qstart, qend = entry.boundaries()
    m_ovl = tree[entry.chrom].overlap(qstart, qend)
    if len(m_ovl) != 1:
        return False
    m_ovl = list(m_ovl)[0]
    end_within = entry.svtype() != truvari.SV.INS

    return truvari.coords_within(qstart, qend, m_ovl.begin, m_ovl.end - 1, end_within)


def entry_within(entry, rstart, rend):
    """
    Extract entry boundaries and type to call `coords_within`
    """
    qstart, qend = entry.boundaries()
    end_within = entry.svtype() != truvari.SV.INS
    return coords_within(qstart, qend, rstart, rend, end_within)


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
