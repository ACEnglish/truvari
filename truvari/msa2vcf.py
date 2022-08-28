"""
Turn an MSA fasta into VCF. Assumes one entry is reference with name >ref_chrom:start-end
"""
import sys
import logging
from io import StringIO
from collections import defaultdict

import pysam

REFIDX = 3
ALTIDX = 4

def build_header(chrom=None):
    """
    Bare minimum vcf header
    """
    ret = '##fileformat=VCFv4.1\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    if chrom is not None:
        ret += f"\n##contig=<ID={chrom}>\n"
    return ret

def msa_to_vars(msa, ref_seq, chrom, start_pos=0):
    """
    Turn MSA into VCF entries and their presence in samples
    returns list of sample names parsed and dictionary of variant : samples containing the variant
    """
    sample_names = set()

    final_vars = defaultdict(list)

    for alt_key in msa.references:
        if alt_key.startswith("ref_"):
            continue
        # gross
        cur_samp_hap = "_".join(alt_key.split('_')[:2])
        sample_names.add(alt_key.split('_')[0])
        alt_seq = msa[alt_key].upper()
        # Gotta assume the first base is a match (little unsafe)
        anchor_base = ref_seq[0]
        if anchor_base == '-':
            logging.error("MSA starts with an indel in %s. Can't make VCF", alt_key)
            sys.exit(1)

        cur_variant = []
        cur_pos = start_pos
        # This is too long. need to have a separate zip method
        for ref_base, alt_base in zip(ref_seq, alt_seq):
            is_ref = ref_base != '-'
            if ref_base == '-':
                ref_base = ""
            if alt_base == '-':
                alt_base = ""

            # nothing to compare
            if not ref_base and not alt_base:
                continue

            if ref_base == alt_base: # No variant
                if cur_variant and is_ref: # back to matching reference
                    # Anchor base correction
                    if len(cur_variant[REFIDX]) == len(cur_variant[ALTIDX]):
                        cur_variant[REFIDX] = cur_variant[REFIDX][1:]
                        cur_variant[ALTIDX] = cur_variant[ALTIDX][1:]
                        cur_variant[1] += 1
                    key = "\t".join([str(_) for _ in cur_variant])
                    # this is a weird edge check
                    # sometimes reference bases aren't aligned
                    if cur_variant[REFIDX] != cur_variant[ALTIDX]:
                        final_vars[key].append(cur_samp_hap)
                    cur_variant = []
            else:
                if not cur_variant:
                    # -1 for the anchor base we're forcing on
                    cur_variant = [chrom, cur_pos - 1, '.', anchor_base + ref_base, anchor_base + alt_base, '.', '.', '.', 'GT']
                else:
                    cur_variant[REFIDX] += ref_base
                    cur_variant[ALTIDX] += alt_base
            if is_ref:
                cur_pos += 1
                anchor_base = ref_base
        # End Zipping
        if cur_variant:
            # Anchor base correction
            if len(cur_variant[REFIDX]) == len(cur_variant[ALTIDX]):
                cur_variant[REFIDX] = cur_variant[REFIDX][1:]
                cur_variant[ALTIDX] = cur_variant[ALTIDX][1:]
                cur_variant[1] += 1
            key = "\t".join([str(_) for _ in cur_variant])
            # this is a weird edge check
            # sometimes reference bases aren't aligned
            if cur_variant[REFIDX] != cur_variant[ALTIDX]:
                final_vars[key].append(cur_samp_hap)
        # End alignment
    sample_names = sorted(list(sample_names))
    return sample_names, final_vars

def make_vcf(variants, sample_names, header):
    """
    Write VCF - building GTs
    """
    out = StringIO()
    out.write(header)
    out.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT\t")
    out.write("\t".join(sample_names) + '\n')

    for var in variants:
        out.write(var)
        for sample in sample_names:
            out.write('\t')
            gt = ["0", "0"]
            if sample + '_1' in variants[var]:
                gt[0] = "1"
            if sample + '_2' in variants[var]:
                gt[1] = "1"
            out.write("/".join(gt))
        out.write('\n')
    out.seek(0)
    return out.read()

def msa2vcf(fn, header=None):
    """
    Parses an MSA and returns its VCF as a string
    Assumes one entry in the MSA has the name ref_${chrom}:${start}-${end}
    """
    msa = pysam.FastaFile(fn)

    # The ref_key identifies reference
    ref_key = [_ for _ in msa.references if _.startswith("ref_")][0] # pylint: disable=not-an-iterable
    ref_seq = msa[ref_key].upper()
    chrom = ref_key.split('_')[1]
    start_pos = int(ref_key.split('_')[2])

    sample_names, variants = msa_to_vars(msa, ref_seq, chrom, start_pos)
    if header is None:
        header = build_header(chrom)

    return make_vcf(variants, sample_names, header)
