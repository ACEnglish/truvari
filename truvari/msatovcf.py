"""
Turn an MSA fasta into VCF. Assumes one entry is reference with name >ref_chrom:start-end
"""
import copy
from io import StringIO
from collections import defaultdict

POSIDX = 1
REFIDX = 3
ALTIDX = 4


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

    # If anchor base is identical, we can move down
    # Stop when there's only one base left or an unmached anchor base
    trim = 0
    while trim < len(ref) - 1 and trim < len(alt) - 1 and ref[trim] == alt[trim]:
        trim += 1
    cur_variant[1] += trim
    cur_variant[REFIDX] = ref[trim:]
    cur_variant[ALTIDX] = alt[trim:]
    if len(cur_variant[REFIDX]) == 1 or len(cur_variant[ALTIDX]) == 1:
        return [var_to_str(cur_variant)]

    # decompose REPL to DEL and INS - easier for truvari to compare
    del_var = copy.copy(cur_variant)
    del_var[ALTIDX] = del_var[ALTIDX][0]
    del_var[REFIDX] = del_var[REFIDX][:-1]
    in_var = copy.copy(cur_variant)
    in_var[REFIDX] = in_var[REFIDX][-1]
    in_var[ALTIDX] = in_var[ALTIDX][1:]
    in_var[POSIDX] += len(cur_variant[REFIDX]) - 1
    return [var_to_str(del_var), var_to_str(in_var)]


def aln_to_vars(chrom, start_pos, ref_seq, alt_seq, anchor_base):
    """
    Zip the bases of an alignment and turn into variants
    """
    cur_variant = []
    cur_pos = start_pos
    # This is too long. need to have a separate zip method
    for ref_base, alt_base in zip(ref_seq, alt_seq):
        is_ref = ref_base != '-'
        if ref_base == '-':
            ref_base = ""
        if alt_base == '-':
            alt_base = ""

        # gap on gap
        if not ref_base and not alt_base:
            continue

        if ref_base == alt_base:  # No variant
            if cur_variant and is_ref:  # back to matching reference
                for variant in decompose_variant(cur_variant):
                    yield variant
                cur_variant = []
        else:
            if not cur_variant:
                # -1 for the anchor base we're forcing on
                cur_variant = [chrom, cur_pos - 1, '.', anchor_base + ref_base,
                               anchor_base + alt_base, '.', '.', '.', 'GT']
            else:
                cur_variant[REFIDX] += ref_base
                cur_variant[ALTIDX] += alt_base
        if is_ref:
            cur_pos += 1
            anchor_base = ref_base
    # End Zipping
    if cur_variant:
        for variant in decompose_variant(cur_variant):
            yield variant

def msa_to_vars(msa, chrom, ref_seq=None, start_pos=0, abs_anchor_base='N'):
    """
    Turn MSA into VCF entries and their presence in samples
    returns list of sample names parsed and dictionary of variant : samples containing the variant
    """
    sample_names = set()
    final_vars = defaultdict(list)
    for alt_key in msa.keys():
        if alt_key.startswith("ref_"):
            continue

        # Trim off the location then haplotype from key sample
        cur_samp_hap = alt_key[:alt_key.rindex('_')]
        sample_names.add(cur_samp_hap[:-2])
        if isinstance(msa[alt_key], tuple):
            ref_seq, alt_seq = msa[alt_key]
            ref_seq = ref_seq.upper()
            alt_seq = alt_seq.upper()
        else:
            alt_seq = msa[alt_key].upper()

        anchor_base = ref_seq[0] if ref_seq[0] != '-' else abs_anchor_base
        for variant in aln_to_vars(chrom, start_pos, ref_seq, alt_seq, anchor_base):
            final_vars[variant].append(cur_samp_hap)
    return sorted(list(sample_names)), final_vars


def make_vcf(variants, sample_names):
    """
    Write VCF lines - building GTs
    """
    out = StringIO()
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


def msa2vcf(msa, anchor_base='N'):
    """
    Parse an MSA dict of {name: alignment, ...} and returns its VCF entries as a string

    Assumes one entry in the MSA has the name `ref_${chrom}:${start}-${end}` which gives VCF entries coordinates
    Provide anchor_base to prevent 'N' from being used as an anchor base
    Returns a string of entries

    Example (for dealing with test coverage not being seen)
        >>> import truvari
        >>> from truvari.phab import fasta_reader
        >>> msa_dir = "repo_utils/test_files/external/fake_mafft/lookup/"
        >>> msa_file = "fm_7bb50c57d657828978076072c80f8a1f.msa"
        >>> seqs = open(msa_dir + msa_file).read()
        >>> fasta = dict(fasta_reader(seqs))
        >>> m_entries_str = truvari.msa2vcf(fasta)
    """
    ref_key = [_ for _ in msa.keys() if _.startswith("ref_")][0]
    chrom, rest = ref_key[len("ref_"):].split(':')
    ref_seq = msa[ref_key].upper() if isinstance(msa[ref_key], str) else None
    start_pos = int(rest.split('-')[0])

    sample_names, variants = msa_to_vars(
        msa, chrom, ref_seq, start_pos, anchor_base)
    return make_vcf(variants, sample_names)
