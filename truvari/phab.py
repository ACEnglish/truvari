"""
Wrapper around MAFFT and other tools to perform an MSA comparison of phased variants
"""
import os
import re
import sys
import shutil
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
    # could maybe allow this to be an existing MSA and then we just use comp and add-to
    parser.add_argument("-c", "--comp", type=str, default=None,
                        help="Comparison vcf to MSA")
    parser.add_argument("-f", "--reference", type=str, required=True,
                        help="Reference")
    parser.add_argument("--buffer", type=int, default=100,
                        help="Number of reference bases before/after region to add to MSA")
    #parser.add_argument("--add-to", action='store_true',
    #                    help="Build the baseMSA independentally of the compMSA, then add the comp")
    parser.add_argument("-o", "--output", default="phab_out",
                        help="Output directory")
    parser.add_argument("--bSamples", type=str, default=None,
                        help="Subset of samples to MSA from base-VCF")
    parser.add_argument("--cSamples", type=str, default=None,
                        help="Subset of samples to MSA from comp-VCF")
    args = parser.parse_args(args)
    return args

def get_reference(fn, chrom, start, end, output):
    """
    Pull a subset of the reference - name it so we can tie back later
    """
    fasta = pysam.FastaFile(fn)
    # Ugh.. not sure about this -1 anymore
    oseq = fasta.fetch(chrom, start - 1, end)
    # no need to make it pretty?
    #oseq = re.sub("(.{60})", "\\1\n", oseq, 0, re.DOTALL)
    with open(output, 'w') as fout:
        fout.write(f">ref_{chrom}_{start}_{end}\n{oseq}\n")

def pull_variants(vcf, region, output, ref, samples=None):
    """
    Given vcf and a region, grab the relevant variants
    Maybe can make fill-from fasta optional.. but it isn't huge overhead (for now)
    """
    if samples is not None:
        samples = "-s " + ",".join(samples)
    else:
        samples = ""
    cmd = f"""bcftools view -c 1 {samples} -r {region} {vcf} \
| bcftools +fill-from-fasta /dev/stdin -- -c REF -f {ref} \
| bgzip > {output} && tabix {output}"""

    ret = truvari.cmd_exe(cmd, pipefail=True)
    if ret.ret_code != 0:
        logging.error("Unable to pull variants from %s", vcf)
        logging.error(ret.stderr)
        sys.exit(1)

def build_consensus(vcf, ref, region, output, samples=None):
    """
    Make the consensus sequence - appends to output
    """
    cmd = """samtools faidx {ref} {region} | \
bcftools consensus -H{hap} --sample {sample} --prefix {sample}_{hap}_ {vcf} >> {output}
"""
    for samp in samples:
        for hap in [1, 2]:
            m_cmd = cmd.format(ref=ref,
                               region=region,
                               hap=hap,
                               sample=samp,
                               vcf=vcf,
                               output=output)
            ret = truvari.cmd_exe(m_cmd, pipefail=True)
            if ret.ret_code != 0:
                logging.error("Unable to make haplotype %d for sample %s", hap, samp)
                logging.error(ret.stderr)
                sys.exit(1)

def run_mafft(seq_fn, output):
    """
    Run mafft
    """
    cmd = f"mafft.bat --retree 2 --maxiterate 0 {seq_fn} > {output}"
    ret = truvari.cmd_exe(cmd)
    if ret.ret_code != 0: # this doesn't (ever?) exit non-zero
        logging.error("Unable to run MAFFT on %s", seq_fn)
        logging.error(ret.stderr)
        sys.exit(1)

def run_phab(args):
    """
    Regularize variants with MSA
    """
    # make buffered region coordinates
    chrom, start, end = re.split(':|-', args.region)
    start = int(start) - args.buffer
    end = int(end) + args.buffer
    buff_region = f"{chrom}:{start}-{end}"

    try:
        os.mkdir(args.output)
    except FileExistsError:
        # Remove this agression
        pass
    # if samples is none, we're building consensus for all the samples
    if args.bSamples is None:
        args.bSamples = list(pysam.VariantFile(args.base).header.samples)
    if args.comp and args.cSamples is None:
        args.cSamples = list(pysam.VariantFile(args.comp).header.samples)
    # need to parse them out here

    sequences = os.path.join(args.output, "haps.fa")
    get_reference(args.reference, chrom, start, end, sequences)

    base_input_vcf = os.path.join(args.output, "base.vcf.gz")
    pull_variants(args.base, args.region, base_input_vcf, args.reference, args.bSamples)
    build_consensus(base_input_vcf, args.reference, buff_region, sequences, args.bSamples)

    # if add-to, we put this elsewhere and revisit for now I'll disable that feature
    # because it assumes we already have an MSA
    if args.comp is not None:
        comp_input_vcf = os.path.join(args.output, "comp.vcf.gz")
        pull_variants(args.comp, args.region, comp_input_vcf, args.reference, args.cSamples)
        build_consensus(args.comp, args.reference, buff_region, sequences, args.cSamples)

    msa_output = os.path.join(args.output, "msa.fa")
    run_mafft(sequences, msa_output)

    output_vcf = os.path.join(args.output, "output.vcf")
    vcf = pysam.VariantFile(args.base)
    n_header = '##fileformat=VCFv4.1\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    n_header += str(vcf.header.contigs[chrom].header_record)

    with open(output_vcf, 'w') as fout:
        fout.write(truvari.msa2vcf(msa_output, n_header))

    truvari.compress_index_vcf(output_vcf)

    # I could do reports here...

def check_requirements():
    """
    ensure external programs are in PATH
    """
    check_fail = False
    for prog in ["bcftools", "bgzip", "tabix", "samtools", "mafft.bat"]:
        if not shutil.which(prog):
            logging.error("Unable to find `%s` in PATH", prog)
            check_fail = True
    return check_fail

def phab_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)
    truvari.setup_logging(True)
    if check_requirements():
        sys.exit(1)
    #if check_params(args):
        # mainly just if files exist
        #pass
    # need to be api-ize
    run_phab(args)
    logging.info("Finished phab")
