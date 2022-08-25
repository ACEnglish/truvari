"""
Collection of methods with event helpers
that compare events, coordinates, or transform vcf entries
"""
import re
import logging

import edlib
import Levenshtein

import truvari


def entry_is_present(entry, sample=None):
    """
    Checks if entry's sample genotype is present and is heterozygous or homozygous (a.k.a. present)

    :param `entry`: entry to check
    :type `entry`: :class:`pysam.VariantRecord`
    :param `sample`: sample name
    :type `sample`: string, optional

    :return: True if variant is present in the sample
    :rtype: bool

    Example
        >>> import truvari
        >>> import pysam
        >>> v = pysam.VariantFile('repo_utils/test_files/input1.vcf.gz')
        >>> truvari.entry_is_present(next(v))
        True
    """
    if sample is None:
        sample = entry.samples.keys()[0]
    return "GT" in entry.samples[sample] and \
           truvari.get_gt(entry.samples[sample]["GT"]) in [
        truvari.GT.HET, truvari.GT.HOM]


def entry_to_key(entry, prefix="", bounds=False):
    """
    Turn a vcf entry into a hashable key string. Use the prefix (base/comp) to indicate the source
    VCF when consolidating multiple files' calls. If bounds: call entry_boundaries for start/stop.

    .. warning::
        If a caller redundantly calls a variant exactly the same. It will not have a unique key

    :param `source`: string to identify vcf from which this entry originates
    :type `source`: string
    :param `entry`: entry to stringify
    :type `entry`: :class:`pysam.VariantRecord`

    :return: hashable string uniquely identifying the variant
    :rtype: string
    """
    if prefix:
        prefix += '.'
    if bounds:
        start, end = entry_boundaries(entry)
        return f"{prefix}{entry.chrom}:{start}-{end}({entry.ref}|{entry.alts[0]})"

    return f"{prefix}{entry.chrom}:{entry.start}-{entry.stop}.{entry.alts[0]}"


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
    return min(sizeA, sizeB) / float(max(sizeA, sizeB)), sizeA - sizeB


def entry_distance(entryA, entryB):
    """
    Calculate the start and end distances of the pair
    """
    astart, aend = entry_boundaries(entryA)
    bstart, bend = entry_boundaries(entryB)
    return astart - bstart, aend - bend


def entry_size_similarity(entryA, entryB):
    """
    Calculate the size similarity percent for two entries

    :param `entryA`: first entry
    :type `entryA`: :class:`pysam.VariantRecord`
    :param `entryB`: second entry
    :type `entryB`: :class:`pysam.VariantRecord`

    :return: size similarity percent and size diff (A - B)
    :rtype: (float, int)

    Example
        >>> import truvari
        >>> import pysam
        >>> v = pysam.VariantFile('repo_utils/test_files/input1.vcf.gz')
        >>> a = next(v)
        >>> b = next(v)
        >>> truvari.entry_size_similarity(a, b)
        (0.07142857142857142, 13)
    """
    sizeA = entry_size(entryA)
    sizeB = entry_size(entryB)
    return sizesim(sizeA, sizeB)


def entry_gt_comp(entryA, entryB, sampleA=None, sampleB=None):
    """
    Compare the genotypes of two entries

    :param `entryA`: first entry
    :type `entryA`: :class:`pysam.VariantRecord`
    :param `entryB`: second entry
    :type `entryB`: :class`pysam.VariantRecord`
    :param `sampleA`: sample of entryA to check
    :type `sampleA`: string, optional
    :param `sampleB`: sample of entryB to check
    :type `sampleB`: string, optional

    :return: True if the genotypes are concordant
    :rtype: bool

    Example
        >>> import truvari
        >>> import pysam
        >>> v = pysam.VariantFile('repo_utils/test_files/input1.vcf.gz')
        >>> a = next(v)
        >>> b = next(v)
        >>> truvari.entry_gt_comp(a, b)
        True
    """
    if not sampleA:
        sampleA = entryA.samples.keys()[0]
    if not sampleB:
        sampleB = entryB.samples.keys()[0]
    return truvari.get_gt(entryA.samples[sampleA]["GT"]) == truvari.get_gt(entryB.samples[sampleB]["GT"])


def create_pos_haplotype(a1, a2, ref, min_len=0):
    """
    Create haplotypes of two allele's regions that are assumed to be overlapping

    :param `a1`: tuple of chrom, start, end, seq
    :type `a1`: tuple
    :param `a2`: tuple of chrom, start, end, seq
    :type `a2`: tuple
    :param `ref`: Reference genome
    :type `ref`: :class:`pysam.FastaFile`
    :param `min_len`: Minimum length of the haplotype sequence to create
    :type `min_len`: int, optional

    :return: allele haplotype sequences created
    :rtupe: tuple (string, string)
    """
    chrom, a1_start, a1_end, a1_seq = a1
    chrom, a2_start, a2_end, a2_seq = a2
    start = min(a1_start, a2_start)
    end = max(a1_end, a2_end)

    hap_len1 = (abs(a1_start - start) + len(a1_seq) + abs(a1_end - end))
    hap_len2 = (abs(a2_start - start) + len(a2_seq) + abs(a2_end - end))
    min_size = min(hap_len1, hap_len2)
    if min_size < min_len:
        start -= (min_len - min_size) // 2
        end += (min_len + min_size) // 2
    # no negative fetch
    start = max(0, start)
    hap1_seq = ref.fetch(chrom, start, a1_start) + \
        a1_seq + ref.fetch(chrom, a1_end, end)
    hap2_seq = ref.fetch(chrom, start, a2_start) + \
        a2_seq + ref.fetch(chrom, a2_end, end)
    return str(hap1_seq), str(hap2_seq)


def entry_create_haplotype(entryA, entryB, ref, use_ref_seq=False, min_len=0):
    """
    Turn two entries into their haplotype sequence for comparison

    :param `entryA`: first entry
    :type `entryA`: :class:`pysam.VariantRecord`
    :param `entryB`: second entry
    :type `entryB`: :class:`pysam.VariantRecord`
    :param `ref`: Reference genome
    :type `ref`: :class:`pysam.FastaFile`
    :param `use_ref_seq`: If True, use the reference genome to get the sequence instead of the vcf entries
    :type `use_ref_seq`: bool, optional
    :param `min_len`: Minimum length of the haplotype sequence to create
    :type `min_len`: int, optional

    :return: allele haplotype sequences created
    :rtype: tuple : (string, string)
    """
    def get_props(entry):
        """
        We compare the longer of the ref/alt sequence to increase comparability
        """
        if use_ref_seq and (entry.alts[0] == "<DEL>" or len(entry.alts[0]) < len(entry.ref)):
            return entry.chrom, entry.start, entry.stop, ref.fetch(entry.chrom, entry.start, entry.stop)
        return entry.chrom, entry.start, entry.stop, entry.alts[0]
    a1 = get_props(entryA)
    a2 = get_props(entryB)
    return create_pos_haplotype(a1, a2, ref, min_len=min_len)

def entry_to_haplotype(entry, ref, start=None, end=None):
    """
    Integrate an entry into the reference and return the sequence.

    :param `entry`: vcf entry
    :type `entry: :class:`pysam.VariantRecord`
    :param `ref`: reference genome
    :type `ref`: :class:`pysam.FastaFile`
    :param `start`: beginning region of reference to use
    :type `start`: int, optional
    :param `end`: beginning region of reference to use
    :type `end`: int, optional

    :return: allele haplotype sequence
    :rtype: str
    """
    chrom = entry.chrom
    a1_start = entry.start
    a1_end = entry.stop
    hap1_seq = ref.fetch(chrom, start, a1_start) + \
        entry.alts[0] + ref.fetch(chrom, a1_end, end)
    return hap1_seq

def entry_pctsim(entryA, entryB, ref, min_len=0, use_lev=True):
    """
    Calculate similarity of two entries' haplotype changes

    :param `entryA`: first entry
    :type `entryA`: :class:`pysam.VariantRecord`
    :param `entryB`: second entry
    :type `entryB`: :class:`pysam.VariantRecord`
    :param `ref`: Reference genome
    :type `ref`: :class:`pysam.FastaFile`
    :param `min_len`: Percent of selected region's range length to buffer
    :type `min_len`: float, optional
    :param `use_lev`: Use levenshtein distance by default. Set to False to use the faster edlib
    :type `use_lev`: bool, optional

    :return: sequence similarity
    :rtype: float
    """
    # Shortcut to save compute - probably unneeded
    if entryA.ref == entryB.ref and entryA.alts[0] == entryB.alts[0]:
        return 1.0

    if entry_variant_type(entryA) == 'INV' and entry_variant_type(entryB) == 'INV':
        allele1 = entryA.alts[0]
        allele2 = entryB.alts[0]
        return seqsim(allele1, allele2, use_lev)

    # Handling of breakends should be here
    allele1, allele2 = entry_create_haplotype(entryA, entryB, ref, min_len=min_len)
    return seqsim(allele1, allele2, use_lev)


def seqsim(allele1, allele2, use_lev=False):
    """
    Calculate similarity of two sequences

    :param `allele1`: first entry
    :type `allele1`: :class:`pysam.VariantRecord`
    :param `allele2`: second entry
    :type `allele2`: :class:`pysam.VariantRecord`
    :param `use_lev`: Use levenshtein distance by default. Set to False to use the faster edlib
    :type `use_lev`: bool, optional

    :return: sequence similarity
    :rtype: float
    """
    allele1 = allele1.upper()
    allele2 = allele2.upper()
    if use_lev:
        return Levenshtein.ratio(allele1, allele2)
    scr = edlib.align(allele1, allele2)
    totlen = len(allele1) + len(allele2)
    return (totlen - scr["editDistance"]) / totlen


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


def entry_variant_type(entry):
    """
    Return entry's SVTYPE

    .. note::
        How svtype is determined:

        - Starts by trying to use INFO/SVTYPE
        - If SVTYPE is unavailable, infer if entry is a insertion or deletion by \
        looking at the REF/ALT sequence size differences
        - If REF/ALT sequences are not available, try to parse the <INS>, <DEL>, \
        etc from the ALT column.
        - Otherwise, assume 'UNK'

    :param `entry`:
    :type `entry`: :class:`pysam.VariantRecord`

    :return: the determined svtype
    :rtype: string
    """
    sv_alt_match = re.compile(r"\<(?P<SVTYPE>.*)\>")

    ret_type = None
    if "SVTYPE" in entry.info:
        ret_type = entry.info["SVTYPE"]
        if isinstance(ret_type, (list, tuple)):
            ret_type = ret_type[0]
        return ret_type

    if not (entry.alts[0].count("<") or entry.alts[0].count(":")):
        # Doesn't have <INS> or BNDs as the alt seq, then we can assume it's sequence resolved..?
        if len(entry.ref) <= len(entry.alts[0]):
            ret_type = "INS"
        elif len(entry.ref) >= len(entry.alts[0]):
            ret_type = "DEL"
        elif len(entry.ref) == len(entry.alts[0]):
            # Is it really?
            ret_type = "COMPLEX"
        return ret_type
    mat = sv_alt_match.match(entry.alts[0])
    if mat is not None:
        return mat.groupdict()["SVTYPE"]
    logging.warning(
        "SVTYPE is undetermined for entry, using 'UNK' - %s", str(entry))
    return "UNK"


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
    a_type = entry_variant_type(entryA)
    b_type = entry_variant_type(entryB)
    if dup_to_ins and a_type == 'DUP':
        a_type = 'INS'
    if dup_to_ins and b_type == 'DUP':
        b_type = 'INS'
    return a_type == b_type


def entry_boundaries(entry, ins_inflate=False):
    """
    Get the start/end of an entry and order (start < end)

    :param `entry`: entry to get bounds
    :type `entry`: :class:`pysam.VariantRecord`
    :param `ins_inflate`: inflate INS boundaries by sv length
    :type `ins_inflate`: bool, optional

    :return: the entry's start/end boundaries
    :rtype: tuple (int, int)
    """
    start = entry.start
    end = entry.stop
    if ins_inflate and entry_variant_type(entry) == 'INS':
        size = entry_size(entry)
        start -= size // 2
        end += size // 2
    return start, end


def entry_size(entry):
    """
    Determine the size of the variant.

    .. note:: How size is determined

        - Starts by trying to use INFO/SVLEN
        - If SVLEN is unavailable and ALT field is an SV (e.g. <INS>, <DEL>, etc), \
        use abs(vcf.start - vcf.end). The INFO/END tag needs to be available, \
        especially for INS.
        - If len of vcf.REF and vcf.ALT[0] are equal, return len of vcf.REF
        - Otherwise, return the size difference of the sequence resolved call using \
        abs(len(vcf.REF) - len(str(vcf.ALT[0])))

    :param `entry`: entry to look at
    :type `entry`: :class:`pysam.VariantRecord`

    :return: the entry's size
    :rtype: int
    """
    if "SVLEN" in entry.info:
        if type(entry.info["SVLEN"]) in [list, tuple]:
            size = abs(entry.info["SVLEN"][0])
        else:
            size = abs(entry.info["SVLEN"])
    elif entry.alts[0].count("<"):
        start, end = entry_boundaries(entry)
        size = end - start
    else:
        r_len = len(entry.ref)
        a_len = len(entry.alts[0])
        if r_len == a_len:
            size = r_len
        else:
            size = abs(r_len - a_len)
    return size


def weighted_score(sim, size, ovl):
    """
    Unite the similarity measures and make a score
    `(2*sim + 1*size + 1*ovl) / 3.0` scaled to 0-100

    :param `sim`:
    :type `sim`: float
        Sequence similarity
    :param `size`:
    :type `size`: float
        Size similarity
    :param `ovl`: Reciprocal overlap
    :type `ovl`: float

    :return: The score
    :rtype: float
    """
    score = (sim + size + ovl) / 3.0 * 100
    #new_score = score / 1.333333 * 100
    return score


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

def entry_reciprocal_overlap(entry1, entry2):
    """
    Calculates reciprocal overlap of two entries

    :param `entry1`: First entry
    :type `entry1`: :class:`pysam.VariantRecord`
    :param `entry2`: Second entry
    :type `entry2`: :class:`pysam.VariantRecord`

    :return: The reciprocal overlap
    :rtype: float

    Example
        >>> import truvari
        >>> import pysam
        >>> v = pysam.VariantFile('repo_utils/test_files/input1.vcf.gz')
        >>> a = next(v)
        >>> b = next(v)
        >>> truvari.entry_reciprocal_overlap(a, b)
        0
    """
    astart, aend = entry_boundaries(entry1, True)
    bstart, bend = entry_boundaries(entry2, True)
    return reciprocal_overlap(astart, aend, bstart, bend)


def entry_is_filtered(entry, values=None):
    """
    Returns if entry should be filtered given the filter values provided.
    If values is None, assume that filter must have PASS or be blank '.'

    :param `entry`: entry to check
    :type `entry`: :class:`pysam.VariantRecord`
    :param `values`: set of filter values for intersection
    :type `values`: set, optional

    :return: True if entry should be filtered
    :rtype: bool

    Example
        >>> import truvari
        >>> import pysam
        >>> v = pysam.VariantFile('repo_utils/test_files/input1.vcf.gz')
        >>> truvari.entry_is_filtered(next(v)) # PASS shouldn't be filtered
        False
        >>> truvari.entry_is_filtered(next(v), set(["lowQ"])) # Call isn't lowQ, so filter
        True
    """
    if values is None:
        return len(entry.filter) != 0 and 'PASS' not in entry.filter
    return len(set(values).intersection(set(entry.filter))) == 0

def unroll_compare(seqA, seqB, p):
    """
    Unroll two sequences and compare.
    See https://gist.github.com/ACEnglish/1e7421c46ee10c71bee4c03982e5df6c for details

    :param `seqA`: sequence
    :param `seqB`: sequence
    :param `p`: positional difference

    :return: sequence similarity of seqA vs seqB after unrolling
    :rtype: float
    """
    f = p % len(seqB)
    uB = seqB[-f:] + seqB[:-f]
    return seqsim(seqA, uB)
