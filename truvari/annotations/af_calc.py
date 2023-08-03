"""
Extra methods for calculating Allele Fequency annotations
Calculations are python implementations of bcftools +fill_tags
(https://samtools.github.io/bcftools/) Results are within 1e-6
difference of bcftools +fill_tags
"""
import numpy as np


def calc_hwe(nref, nalt, nhet):
    """
    Calculate Hardy Weinberg equilibrium and excess heterozygosity

    :param `nref`: Number of reference alleles
    :type `nref`: int
    :param `nalt`: Number of alternate alleles
    :type `nalt`: int
    :param `nhet`: Number of heterozygous sites
    :type `nhet`: int

    :return: The HWE and ExcHet calculated
    :rtype: tuple (float, float)
    """
    # void calc_hwe(args_t *args, int nref, int nalt, int nhet, float *p_hwe, float *p_exc_het)
    ngt = (nref + nalt) // 2  # also, just nsamples
    nrare = min(nref, nalt)
    # We got an hts_expand array for het probabilities
    probs = np.zeros(nrare + 1, dtype=np.dtype('d'))

    # start at midpoint
    mid = nrare * (nref + nalt - nrare) // (nref + nalt)

    # check to ensure that midpoint and rare alleles have same parity
    if (nrare & 1) ^ (mid & 1):
        mid += 1

    het = mid
    hom_r = (nrare - mid) // 2
    hom_c = ngt - het - hom_r
    probs[mid] = 1
    my_sum = 1

    while het > 1:
        probs[het - 2] = probs[het] * het * \
            (het - 1.0) / (4.0 * (hom_r + 1.0) * (hom_c + 1.0))
        my_sum += probs[het - 2]

        # 2 fewer heterozygotes for next iteration -> add rare, one common homozygote
        hom_r += 1
        hom_c += 1
        # loop
        het -= 2

    het = mid
    hom_r = (nrare - mid) // 2
    hom_c = ngt - het - hom_r
    while het <= nrare - 2:
        probs[het + 2] = probs[het] * 4.0 * hom_r * \
            hom_c / ((het + 2.0) * (het + 1.0))
        my_sum += probs[het + 2]

        # add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
        hom_r -= 1
        hom_c -= 1
        # loop
        het += 2

    probs = probs / my_sum

    p_exc_het = probs[nhet:].sum()
    p_hwe = min(probs[probs > probs[nhet]].sum(), 1)
    return p_exc_het, 1 - p_hwe


def calc_af(gts):
    """
    Calculate allele annotations for a list of genotypes

    :param `gts`: Genotype tuples
    :type `gts`: list

    :return: | Dictonary of
             | AF - allele frequency
             | MAF - minor allele frequency
             | ExcHet - excess heterozygosity
             | HWE - hardy weinberg equilibrium
             | AC - allele count for GT 1
             | MAC - minor allele count
             | AN - number of called alleles
             | N_HEMI - Number of partial genotypes (length of 1 or a single missing)
             | N_HOMREF - homozygous reference GT count
             | N_HET - heterozygous GT count
             | N_HOMALT - homozygous alternate GT count
             | N_MISS - Number of missing genotypes (all .)
    :rtype: dict
    """
    ret = {"AF": 0, "MAF": 0, "ExcHet": 0, "HWE": 0, "MAC": 0, "AC": [0, 0], "AN": 0,
           "N_HEMI": 0, "N_HOMREF": 0, "N_HET": 0, "N_HOMALT": 0, "N_MISS": 0}
    for g in gts:
        if len(g) > 2:
            continue
        for j in g:
            if j is not None:
                ret["AN"] += 1
                ret["AC"][j] += 1
        if len(g) == 2:
            if g[0] == g[1]:
                if g[0] == 0:
                    ret["N_HOMREF"] += 1
                elif g[0] is None:
                    ret["N_MISS"] += 1
                else:
                    ret["N_HOMALT"] += 1
            else:
                if g[0] is None or g[1] is None:
                    ret["N_HEMI"] += 1
                else:
                    ret["N_HET"] += 1
        else:
            if g[0] is None:
                ret["N_MISS"] += 1
            else:
                ret["N_HEMI"] += 1

    if ret["AN"] == 0:
        return ret
    ret["AF"] = ret["AC"][1] / ret["AN"]
    ret["MAC"] = min(ret["AC"])
    ret["MAF"] = ret["MAC"] / ret["AN"]

    ret["ExcHet"], ret["HWE"] = calc_hwe(
        ret["AC"][0], ret["AC"][1], ret["N_HET"])
    return ret


def allele_freq_annos(entry, samples=None):
    """
    Calculate allele annotations for a VCF Entry

    :param `entry`: Entry with samples to parse
    :type `entry`: :class:`pysam.VariantRecord`
    :param `samples`: Subset of samples from the entry over which to calculate annos
    :type `samples`: list of strings, optional

    :return: | Dictonary of
             | AF - allele frequency
             | MAF - minor allele frequency
             | ExcHet - excess heterozygosity
             | HWE - hardy weinberg equilibrium
             | AC - allele count
             | MAC - minor allele count
             | AN - number of called alleles
    :rtype: dict

    Example
        >>> import truvari
        >>> import pysam
        >>> v = pysam.VariantFile('repo_utils/test_files/variants/multi.vcf.gz')
        >>> truvari.allele_freq_annos(next(v))
        {'AF': 0.5, 'MAF': 0.5, 'ExcHet': 1.0, 'HWE': 1.0, 'MAC': 1, 'AC': [1, 1], 'AN': 2, 'N_HEMI': 0, 'N_HOMREF': 0, 'N_HET': 1, 'N_HOMALT': 0, 'N_MISS': 2}
    """
    if samples is None:
        samples = list(entry.samples.keys())
    gts = [entry.samples[_]["GT"] for _ in samples]
    return calc_af(gts)
