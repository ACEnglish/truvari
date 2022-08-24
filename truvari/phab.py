"""
Wrapper around MAFFT and other tools to perform an MSA comparison of phased variants
"""
import re
import logging
import argparse

import pysam
import truvari

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="phab", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-r", "--region", type=str, required=True,
                        help="chrom:start-end of region to process")
    parser.add_argument("-b", "--base", type=str, required=True,
                        help="Baseline vcf to MSA")
    parser.add_argument("-c", "--comp", type=str, default=None,
                        help="Comparison vcf to MSA")
    parser.add_argument("-f", "--reference", type=str, required=True,
                        help="Reference")
    parser.add_argument("--buffer", type=int, default=100,
                        help="Number of reference bases before/after region to add to MSA")
    parser.add_argument("--add-to", action='store_true',
                        help="Build the baseMSA independentally of the compMSA, then add the comp")
    parser.add_argument("-o", "--output", default="phab_out",
                        help="Output directory")
    parser.add_argument("-s", "--samples", type=str, default=None,
                        help="Subset of samples to MSA from base-VCF")
    args = parser.parse_args(args)
    return args

def pull_variants(vcf, output, samples=None):
    """
    given vcf and a region, grab only the relevant variants
    bcftools view -c 1 -r $region $vcf \
    | bcftools +fill-from-fasta /dev/stdin -- -c REF -f $ref \
    | bgzip > msa_${region}/variants.vcf.gz
    """


def get_reference():
    """
    Pul the reference sequence of the region.
    rename contig so we can use it downstream
    # gotta be careful here? I think if I use pysam I'm going to undo my 0/1 based corrections
    """
    pass

def build_consensus():
    """
    Make the consensus sequence
    samtools faidx $ref $expand_region \
        | bcftools consensus -H1 --sample $i --prefix $i_[1|2]_ msa_${region}/variants.vcf.gz \
        | python $DIR/fa_rename.py ${i}_1 >> msa_${region}/haps.fa
    """
    pass

def run_mafft()
    """
    /users/u233287/scratch/misc_software/mafft-linux64/mafft.bat --retree 2 --maxiterate 0 msa_${region}/haps.fa > msa_${region}/aln_results.txt
    """
    pass

def msa2vcf():
    """
    python $DIR/msa2vcf.py msa_${region}/aln_results.txt \
        | vcf-sort \
        | bcftools +fill-tags | bgzip > msa_${region}/result.vcf.gz
    tabix msa_${region}/result.vcf.gz
    """
    pass

def reports():
    pass


def get_reference(fn, chrom, start, end):
    fasta = pysam.FastaFile(fn)
    # Ugh.. not sure about this -1 anymore
    oseq = fasta.fetch(chrom, start - 1, end)
    # no need to make it pretty?
    #oseq = re.sub("(.{60})", "\\1\n", oseq, 0, re.DOTALL)
    print(f">ref_{chrom}_{start}_{end}\n{oseq}")

def run_phab(args):
    
    # make buffered region coordinates
    chrom, start, end = re.split(':|-', args.region)
    start = int(start) - args.buffer
    end = int(end) + args.buffer
    
    pull_variants(args)
    get_reference()
    build_consensus
    run_mafft
    truvari.msa2vcf #

def phab_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)
    if check_requirements(args):
        # samtools bcftools tabix vcf-sort mafft.bat bgzip
        pass
    if check_params(args):
        # mainly just if files exist
        pass
    
    run_phab(args)
