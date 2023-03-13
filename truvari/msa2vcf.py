"""
Turn an MSA fasta into VCF. Assumes one entry is reference with name >ref_chrom:start-end
"""
import copy
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

def decompose_variant(cur_variant):
    """
    Left trim and decompose (repl -> indel) variant
    returns a list of new variants
    """
    def var_to_str(v):
        return "\t".join([str(_) for _ in v])
    ref = cur_variant[REFIDX]
    alt = cur_variant[ALTIDX]
    if ref == alt:
        return []
    # If base after anchor base is identical, we have a new anchor base
    # Stop when there's only one base left or unmached bases
    trim = 0
    while trim < len(ref) - 1 and trim < len(alt) - 1 and ref[trim] == alt[trim]:
        trim += 1
    cur_variant[1] += trim
    cur_variant[REFIDX] = ref[trim:]
    cur_variant[ALTIDX] = alt[trim:]
    if len(cur_variant[REFIDX]) == 1 or len(cur_variant[ALTIDX]) == 1:
        return [var_to_str(cur_variant)]

    # decompose REPL to INS and DEL - easier for truvari to compare
    in_var = copy.copy(cur_variant)
    in_var[REFIDX] = in_var[REFIDX][0]
    del_var = copy.copy(cur_variant)
    del_var[ALTIDX] = del_var[ALTIDX][0]
    return [var_to_str(in_var), var_to_str(del_var)]

def msa_to_vars(msa, ref_seq, chrom, start_pos=0, abs_anchor_base='N'):
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
        cur_samp_hap = "_".join(alt_key.split('_')[:-1])
        sample_names.add("_".join(alt_key.split('_')[:-2]))
        alt_seq = msa[alt_key].upper()

        anchor_base = ref_seq[0] if ref_seq[0] != '-' else abs_anchor_base

        cur_variant = []
        cur_pos = start_pos # still have a problem here with anchor base.
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
                    for key in decompose_variant(cur_variant):
                        final_vars[key].append(cur_samp_hap)
                    cur_variant = []
            else:
                if not cur_variant:
                    # -1 for the anchor base we're forcing on
                    cur_variant = [chrom, cur_pos - 1, '.', anchor_base + ref_base, \
                                   anchor_base + alt_base, '.', '.', '.', 'GT']
                else:
                    cur_variant[REFIDX] += ref_base
                    cur_variant[ALTIDX] += alt_base
            if is_ref:
                cur_pos += 1
                anchor_base = ref_base
        # End Zipping
        if cur_variant:
            for key in decompose_variant(cur_variant):
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

def msa2vcf(fn, header=None, anchor_base=None):
    """
    Parse an MSA file and returns its VCF as a string

    Assumes one entry in the MSA has the name `ref_${chrom}:${start}-${end}` which gives VCF entries coordinates

    VCF can be given a header or a minimal one is created with only file format and a single contig line
    Provide anchor_base to prevent 'N' from being used as an anchor base
    Returns a string
    """
    msa = pysam.FastaFile(fn)

    # The ref_key identifies reference
    ref_key = [_ for _ in msa.references if _.startswith("ref_")][0] # pylint: disable=not-an-iterable
    ref_seq = msa[ref_key].upper()
    _, chrom, start_pos, _ = ref_key.split('_')
    start_pos = int(start_pos)

    sample_names, variants = msa_to_vars(msa, ref_seq, chrom, start_pos, anchor_base)
    if header is None:
        header = build_header(chrom)

    return make_vcf(variants, sample_names, header)
