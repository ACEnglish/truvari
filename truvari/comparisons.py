"""
Collection of methods with event helpers
that compare events, coordinates, or transform vcf entries
"""
# pylint: disable=unused-variable
import re
import logging
from functools import cmp_to_key

import edlib
import Levenshtein

import truvari


def entry_is_present(entry, sample=None):
    """
    Checks if entry's sample genotype is present and is heterozygous or homozygous (a.k.a. present)

    Parameters
    ----------
    entry : pysam.VariantRecord
        entry to check
    sample : string, optional
        sample name

    Returns
    -------
    boolean
        True if variant is present in the sample
    """
    if sample is None:
        sample = entry.samples.keys()[0]
    return "GT" in entry.samples[sample] and \
           truvari.get_gt(entry.samples[sample]["GT"]) in [
        truvari.GT.HET, truvari.GT.HOM]


def entry_to_key(source, entry):
    """
    Turn a vcf entry into a hashable key using the 'source' (base/comp) to separate the two
    helpful for not re-using variants
    BUG: if a caller redundantly calls a variant exactly the same. It will be collapsed here

    Parameters
    ----------
    source : string
        string to identify vcf from which this entry originates
    entry : pysam.VariantRecord
        entry to stringify

    Returns
    -------
    string : hashable string uniquely identifying the variant
    """
    start, end = entry_boundaries(entry)
    return "%s.%s:%d-%d(%s|%s)" % (source, entry.chrom, start, end, entry.ref, entry.alts[0])  # pylint: disable=consider-using-f-string


def sizesim(sizeA, sizeB):
    """
    Calculate the size similarity percent for two sizes

    Parameters
    ----------
    sizeA : int
        first size
    sizeB : int
        second size

    Returns
    -------
    float : size similarity percent
    """
    return min(sizeA, sizeB) / float(max(sizeA, sizeB)), sizeA - sizeB


def entry_size_similarity(entryA, entryB):
    """
    Calculate the size similarity percent for two entries

    Parameters
    ----------
    entryA : pysam.VariantRecord
        first entry
    entryB : pysam.VariantRecord
        second entry

    Returns
    -------
    float : size similarity percent
    """
    sizeA = entry_size(entryA)
    sizeB = entry_size(entryB)
    return sizesim(sizeA, sizeB)


def entry_gt_comp(entryA, entryB, sampleA=None, sampleB=None):
    """
    Compare the genotypes of two entries

    Parameters
    ----------
    entryA : pysam.VariantRecord
        first entry
    entryB : pysam.VariantRecord
        second entry
    sampleA : string, optional
        sample of entryA to check
    sampleB: string, optional
        sample of entryB to check

    Returns
    -------
    boolean
    """
    if not sampleA:
        sampleA = entryA.samples.keys()[0]
    if not sampleB:
        sampleB = entryB.samples.keys()[0]
    return entryA.samples[sampleA]["GT"] == entryB.samples[sampleB]["GT"]


def create_pos_haplotype(chrom, a1_start, a1_end, a1_seq, a2_start, a2_end, a2_seq, ref, buf_len=0):
    """
    Create haplotypes of two allele's regions that are assumed to be overlapping

    Parameters
    ----------
    chrom : string
        Chromosome where alleles are located
    a1_start : integer
        0-based position of first allele's start
    a1_end : integer
        0-based position of first allele's  end
    a1_seq : string
        First allele's sequene
    a2_start : integer
        0-based position of second allele's start
    a2_end : integer
        0-based position of second allele's  end
    a2_seq : string
        Second allele's sequene
    ref : pysam.FastaFile
        Reference genome
    buf_len : float, optional
        Percent of selected region's range length to buffer

    Returns
    -------
    tuple : (string, string)
        allele haplotype sequences created
    """
    start = min(a1_start, a2_start)
    end = max(a1_end, a2_end)
    buff = int((end - start) * buf_len)
    start -= buff
    end += buff
    hap1_seq = ref.fetch(chrom, start, a1_start) + \
        a1_seq + ref.fetch(chrom, a1_end, end)
    hap2_seq = ref.fetch(chrom, start, a2_start) + \
        a2_seq + ref.fetch(chrom, a2_end, end)
    return str(hap1_seq), str(hap2_seq)


def entry_create_haplotype(entryA, entryB, ref, use_ref_seq=False, buf_len=0):
    """
    Turn two entries into their haplotype sequence for comparison

    Parameters
    ----------
    entryA : pysam.VariantRecord
        first entry
    entryB : pysam.VariantRecord
        second entry
    ref : pysam.FastaFile
        Reference genome
    use_ref_seq : boolean, optional
        If True, use the reference genome to get the sequence instead of the vcf entries
    buf_len : float, optional
        Percent of selected region's range length to buffer

    Returns
    -------
    tuple : (string, string)
        allele haplotype sequences created
    """
    def get_props(entry):
        """
        We compare the longer of the ref/alt sequence to increase comparability
        """
        if use_ref_seq and (entry.alts[0] == "<DEL>" or len(entry.alts[0]) < len(entry.ref)):
            return entry.chrom, entry.start, entry.stop, ref.fetch(entry.chrom, entry.start, entry.stop)
        return entry.chrom, entry.start, entry.stop, entry.alts[0]

    a1_chrom, a1_start, a1_end, a1_seq = get_props(entryA)
    a2_chrom, a2_start, a2_end, a2_seq = get_props(entryB)
    return create_pos_haplotype(a1_chrom, a1_start, a1_end, a1_seq, a2_start, a2_end, a2_seq, ref, buf_len=buf_len)


def entry_pctsim(entryA, entryB, ref, buf_len=0, use_lev=True):
    """
    Calculate similarity of two sequences' haplotype changes

    Parameters
    ----------
    entryA : pysam.VariantRecord
        first entry
    entryB : pysam.VariantRecord
        second entry
    ref : pysam.FastaFile
        Reference genome
    buf_len : float, optional
        Percent of selected region's range length to buffer
    use_lev : boolean, optional
        Use levenshtein distance by default. Set to False to use the faster edlib

    Returns
    -------
    float : sequence similarity
    """
    # Shortcut to save compute - probably unneeded
    if entryA.ref == entryB.ref and entryA.alts[0] == entryB.alts[0]:
        return 1.0
    # Handling of breakends should be here
    try:
        allele1, allele2 = entry_create_haplotype(
            entryA, entryB, ref, buf_len=buf_len)
    except Exception as e:  # pylint: disable=broad-except
        logging.critical('Unable to compare sequence similarity\n%s\n%s\n%s', str(
            entryA), str(entryB), str(e))
        return 0
    if use_lev:
        return Levenshtein.ratio(allele1, allele2)
    scr = edlib.align(allele1, allele2)
    totlen = len(allele1) + len(allele2)
    return (totlen - scr["editDistance"]) / totlen


def overlaps(s1, e1, s2, e2):
    """
    Check if two start/end ranges have overlap

    Parameters
    ----------
    s1 : int
        range 1 start
    e1 : int
        range 1 end
    s2 : int
        range 2 start
    e2 : int
        range 2 end

    Returns
    -------
    boolean
    """
    s_cand = max(s1, s2)
    e_cand = min(e1, e2)
    return s_cand < e_cand


def entry_variant_type(entry):
    """
    How svtype is determined:
    - Starts by trying to use INFO/SVTYPE
    - If SVTYPE is unavailable, infer if entry is a insertion or deletion by
      looking at the REF/ALT sequence size differences
    - If REF/ALT sequences are not available, try to parse the <INS>, <DEL>,
      etc from the ALT column.
    - Otherwise, assume 'UNK'

    Parameters
    ----------
    entry : pysam.VariantRecord

    Returns
    -------
    string : the determined svtype
    """
    sv_alt_match = re.compile(r"\<(?P<SVTYPE>.*)\>")

    ret_type = None
    if "SVTYPE" in entry.info:
        ret_type = entry.info["SVTYPE"]
        if isinstance(ret_type, list):
            logging.warning("SVTYPE is list for entry %s", str(entry))
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


def entry_same_variant_type(entryA, entryB):
    """
    returns if entryA svtype == entryB svtype

    Parameters
    ----------
    entryA : pysam.VariantRecord
        first entry
    entryB : pysam.VariantRecord
        second entry

    Returns
    -------
    boolean
    """
    a_type = entry_variant_type(entryA)
    b_type = entry_variant_type(entryB)
    return a_type == b_type


def fetch_coords(lookup, entry, dist=0):
    """
    Get the minimum/maximum fetch coordinates to find all variants within dist of variant

    Parameters
    ----------
    lookup : dict
        genome tree build defaultdict with interevaltrees
    entry : pysam.VariantRecord
        entry to build coords from
    dist : integer
        distance buffer to add/subtract from the coords

    Returns
    -------
    tuple : (int, int)
        the minimum/maximum fetch coordinates for the entry
    """
    start, end = entry_boundaries(entry)
    start -= dist
    end += dist
    # Membership queries are fastest O(1)
    if not lookup[entry.chrom].overlaps(start, end):
        return None, None

    cand_intervals = lookup[entry.chrom].overlap(start, end)
    s_ret = min(
        [x.data for x in cand_intervals if overlaps(start, end, x[0], x[1])])
    e_ret = max(
        [x.data for x in cand_intervals if overlaps(start, end, x[0], x[1])])
    return s_ret, e_ret


def entry_boundaries(entry):
    """
    Get the start/end of an entry and order (start < end)

    Parameters
    ----------
    entry : pysam.VariantRecord
        entry to get bounds

    Returns
    -------
    tuple : (int, int)
        the entry's start/end boundaries
    """
    start = entry.start
    end = entry.stop
    return start, end


def entry_size(entry):
    """
    Determine the size of the variant.

    How size is determined:
    - Starts by trying to use INFO/SVLEN
    - If SVLEN is unavailable and ALT field is an SV (e.g. <INS>, <DEL>, etc),
      use abs(vcf.start - vcf.end). The INFO/END tag needs to be available,
      especially for INS.
    - Otherwise, return the size difference of the sequence resolved call using
      abs(len(vcf.REF) - len(str(vcf.ALT[0])))

    Parameters
    ----------
    entry : pysam.VariantRecord
        entry to look at

    Returns
    -------
    int :
        the entry's size
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
        size = abs(len(entry.ref) - len(entry.alts[0]))
    return size


def weighted_score(sim, size, ovl):
    """
    Unite the similarity measures and make a score
    return (2*sim + 1*size + 1*ovl) / 3.0 scaled to 0-100

    Parameters
    ----------
    sim : float
        Sequence similarity
    size : float
        Size similarity
    ovl : float
        Reciprocal overlap

    Returns
    -------
    float :
        The score
    """
    score = (2 * sim + 1 * size + 1 * ovl) / 3.0
    new_score = int(score / 1.333333 * 100)
    return new_score


def reciprocal_overlap(astart, aend, bstart, bend):
    """
    Calculates a reciprocal overlap for two ranges

    Parameters
    ----------
    astart : int
        First range's start position
    aend : ing
        First range's end position
    bstart : int
        Second range's start position
    bend : int
        Second range's end position

    Returns
    -------
    float :
        The reciprocal overlap
    """
    ovl_start = max(astart, bstart)
    ovl_end = min(aend, bend)
    if ovl_start < ovl_end:  # Otherwise, they're not overlapping
        ovl_pct = float(ovl_end - ovl_start) / \
            max(aend - astart, bend - bstart)
    else:
        ovl_pct = 0
    return ovl_pct


def entry_reciprocal_overlap(entry1, entry2):
    """
    Calculates a reciprocal overlap rule for two entries

    Parameters
    ----------
    entry1 : pysam.VariantRecord
        First entry
    entry2 : pysam.VariantRecord
        Second entry

    Returns
    -------
    float :
        The reciprocal overlap
    """
    astart, aend = entry_boundaries(entry1)
    bstart, bend = entry_boundaries(entry2)
    return reciprocal_overlap(astart, aend, bstart, bend)


def filter_value(entry, values=None):
    """
    Returns if entry should be filtered given the filter values provided.
    If values is None, assume that filter must have PASS or be blank '.'

    Parameters
    ----------
    entry : pysam.VariantRecord
        entry to check

    Returns
    -------
    bool :
        True if entry should be filtered
    """
    if values is None:
        return len(entry.filter) != 0 and 'PASS' not in entry.filter
    return values.intersection(entry.filter)


def match_sorter(candidates):
    """
    Sort a list of MATCHRESULT tuples inplace.

    Parameters
    ----------
    candidates : list of MATCHRESULT named tuples
    """
    if len(candidates) == 0:
        return
    entry_idx = len(candidates[0]) - 1

    def sort_cmp(mat1, mat2):
        """
        Sort by attributes and then deterministically by hash(str(VariantRecord))
        """
        for i in range(entry_idx):
            if mat1[i] != mat2[i]:
                return mat1[i] - mat2[i]
        return hash(str(mat1[entry_idx])) - hash(str(mat2[entry_idx]))

    candidates.sort(reverse=True, key=cmp_to_key(sort_cmp))


def copy_entry(entry, header):
    """
    Make a pysam.VariantRecord editable

    Parameters
    ----------
    entry : pysam.VariantRecord
        entry to make editable
    header : pysam.VariantHeader
        header of output vcf

    Returns
    -------
    pysam.VariantRecord : editable record
    """
    try:
        ret = header.new_record(contig=entry.chrom, start=entry.start, stop=entry.stop,
                                alleles=entry.alleles, id=entry.id, qual=entry.qual, filter=entry.filter,
                                info=entry.info)
    except TypeError as e:
        new_entry_info = dict(entry.info)
        for key, value in new_entry_info.items():
            if isinstance(value, tuple):
                if header.info[key].type != "String":
                    logging.error("Entry is not copyable by pysam. INFO %s has Number=%s and Type=%s",
                                  key, header.info[key].number, header.info[key].type)
                    logging.error(
                        "Number should be changed to '.' or Type to 'String'")
                    logging.error("Check VCF header (%s)", str(entry))
                    raise e
                new_entry_info[key] = ",".join(value)
        ret = header.new_record(contig=entry.chrom, start=entry.start, stop=entry.stop,
                                alleles=entry.alleles, id=entry.id, qual=entry.qual, filter=entry.filter,
                                info=new_entry_info)
    for sample in entry.samples:
        for k, v in entry.samples[sample].items():
            # this will be a problem for pVCFs with differing Number=./A/G and set on input as (None,).. maybe
            try:
                ret.samples[sample][k] = v
            except TypeError:
                pass

    return ret
