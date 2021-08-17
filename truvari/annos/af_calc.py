"""
Extra methods for calculating Allele Fequency annotations
Calculations are python implementations of bcftools +fill_tags
(https://samtools.github.io/bcftools/) Results are within 1e-6
difference of bcftools +fill_tags
"""
from collections import Counter

import numpy as np

def calc_hwe(nref, nalt, nhet):
    """
    Returns tuple of the HWE and the ExeHet
    """
    #void calc_hwe(args_t *args, int nref, int nalt, int nhet, float *p_hwe, float *p_exc_het)
    ngt = (nref + nalt) // 2 # also, just nsamples
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
        probs[het - 2] = probs[het] * het * (het - 1.0) / (4.0 * (hom_r + 1.0) * (hom_c + 1.0))
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
        probs[het + 2] = probs[het] * 4.0 * hom_r * hom_c / ((het + 2.0) * (het + 1.0))
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

def allele_freq_annos(entry, samples=None):
    """
    Given a pysam entry, return dictonary with keys
        AF, MAF, ExcHet, HWE, AC, MAC
    samples is a list of samples to subset from the entry to calculate the allele frequency
    """
    if samples is None:
        samples = list(entry.samples.keys())

    n_samps = 0
    n_het = 0
    cnt = Counter() # 0 or 1 allele counts
    for samp in samples:
        dat = entry.samples[samp]
        if None in dat["GT"]:
            continue
        n_samps += 1
        for j in dat["GT"]:
            cnt[j] += 1
        n_het += 1 if dat["GT"][0] != dat["GT"][1] else 0

    if n_samps == 0:
        return {"AF":0, "MAF":0, "ExcHet":0, "HWE":0, "MAC":0, "AC": 0}

    af = cnt[1] / (n_samps * 2)
    srt = [(v, k) for k, v in sorted(cnt.items(), key=lambda item: item[1])]
    ac = [cnt[_] for _ in [0, 1]]
    mac = srt[0][0]
    maf = 1 - (srt[-1][0] / (n_samps * 2))

    p_exc_het, p_hwe = calc_hwe(cnt[0], cnt[1], n_het)
    return {"AF":af, "MAF":maf, "ExcHet":p_exc_het, "HWE":p_hwe, "MAC":mac, "AC": ac}
