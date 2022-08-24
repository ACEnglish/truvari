"""
Turn an MSA fasta into VCF. Assumes one entry is reference with name >ref_chrom:start-end
"""
import sys
import argparse
import pysam
from collections import defaultdict
import logging
import truvari

REFIDX = 3
ALTIDX = 4

def build_header():
    """
    Pull header contig information from existing VCF if provided, otherwise we'll fake it
    """
    print("##fileformat=VCFv4.1")
    print("##contig=<ID=chr1,length=248956422>")
    print("##contig=<ID=chr2,length=242193529>")
    print("##contig=<ID=chr3,length=198295559>")
    print("##contig=<ID=chr4,length=190214555>")
    print("##contig=<ID=chr5,length=181538259>")
    print("##contig=<ID=chr6,length=170805979>")
    print("##contig=<ID=chr7,length=159345973>")
    print("##contig=<ID=chr8,length=145138636>")
    print("##contig=<ID=chr9,length=138394717>")
    print("##contig=<ID=chr10,length=133797422>")
    print("##contig=<ID=chr11,length=135086622>")
    print("##contig=<ID=chr12,length=133275309>")
    print("##contig=<ID=chr13,length=114364328>")
    print("##contig=<ID=chr14,length=107043718>")
    print("##contig=<ID=chr15,length=101991189>")
    print("##contig=<ID=chr16,length=90338345>")
    print("##contig=<ID=chr17,length=83257441>")
    print("##contig=<ID=chr18,length=80373285>")
    print("##contig=<ID=chr19,length=58617616>")
    print("##contig=<ID=chr20,length=64444167>")
    print("##contig=<ID=chr21,length=46709983>")
    print("##contig=<ID=chr22,length=50818468>")
    print("##contig=<ID=chrX,length=156040895>")
    print("##contig=<ID=chrY,length=57227415>")
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')

def msa_to_vars(msa, startpos=0):
    """
    build dictionary of variant : samples containing the variant
    """
    defaultdict(list)
    for alt_key in msa.references:
        if alt_key.startswith("ref_"):
            continue
        cur_samps = alt_key.split(';')
        sample_names.update([_.split('_')[0] for _ in cur_samps])
        alt_seq = msa[alt_key].upper()
        # Gotta assume the first base is a match (little unsafe)
        anchor_base = ref_seq[0]
        cur_variant = None
        cur_pos = start_pos
        for pos, bases in enumerate(zip(ref_seq, alt_seq)):
            ref_base, alt_base = bases
            is_ref = ref_base != '-'
            if ref_base == '-':
                ref_base = ""
            if alt_base == '-':
                alt_base = ""

            # nothing to compare
            if not ref_base and not alt_base:
                continue

            if ref_base == alt_base: # No variant
                if cur_variant is not None and is_ref: # back to matching reference
                    # Anchor base correction
                    if len(cur_variant[REFIDX]) == len(cur_variant[ALTIDX]):
                        cur_variant[REFIDX] = cur_variant[REFIDX][1:]
                        cur_variant[ALTIDX] = cur_variant[ALTIDX][1:]
                        cur_variant[1] += 1
                    key = "\t".join([str(_) for _ in cur_variant])
                    # this is a weird edge check
                    # sometimes reference bases aren't aligned
                    if cur_variant[REFIDX] != cur_variant[ALTIDX]:
                        final_vars[key].extend(cur_samps)
                    cur_variant = None
            else:
                if cur_variant is None:
                    # -1 for the anchor base we're forcing on
                    cur_variant = [chrom, cur_pos - 1, '.', anchor_base + ref_base, anchor_base + alt_base, '.', '.', '.', 'GT']
                else:
                    cur_variant[REFIDX] += ref_base
                    cur_variant[ALTIDX] += alt_base
            if is_ref:
                cur_pos += 1
                anchor_base = ref_base
        if cur_variant is not None:
            # Anchor base correction
            if len(cur_variant[REFIDX]) == len(cur_variant[ALTIDX]):
                cur_variant[REFIDX] = cur_variant[REFIDX][1:]
                cur_variant[ALTIDX] = cur_variant[ALTIDX][1:]
                cur_variant[1] += 1
            key = "\t".join([str(_) for _ in cur_variant])
            # this is a weird edge check
            # sometimes reference bases aren't aligned
            if cur_variant[REFIDX] != cur_variant[ALTIDX]:
                final_vars[key].extend(cur_samps)
            cur_variant = None

def write_vcf(variants, sample_names):
    """
    Write VCF - building GTs
    """
    sys.stdout.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT\t")
    sys.stdout.write("\t".join(sample_names) + '\n')

    for var in variants:
        sys.stdout.write(var)
        for sample in sample_names:
            sys.stdout.write('\t')
            gt = ["0", "0"]
            if sample + '_1' in final_vars[var]:
                gt[0] = "1"
            if sample + '_2' in final_vars[var]:
                gt[1] = "1"
            sys.stdout.write("/".join(gt))
        sys.stdout.write('\n')

def msa2vcf(fn, ref_prefix="ref_", header=None):
    
    pass

def msa2vcf_main():
    
    msa = pysam.FastaFile(sys.argv[1])

    # The ref_key identifies reference
    ref_key = [_ for _ in msa.references if _.startswith("ref_")][0]
    ref_seq = msa[ref_key].upper()
    chrom = ref_key.split('_')[1]
    start_pos = int(ref_key.split('_')[2])

    # We merge variants that are identical
    final_vars = msa_to_vars() 
    
    sample_names = set()

    sample_names = sorted(list(sample_names))
    write_vcf(final_vars, sample_names)
    logging.info("Finished")
