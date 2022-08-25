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
                        help="Bed filename or comma-separated list of chrom:start-end regions to process")
    # could maybe allow this to be an existing MSA and then we just use comp and add-to
    parser.add_argument("-b", "--base", type=str, required=True,
                        help="Baseline vcf to MSA")
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
    parser.add_argument("-m", "--mafft-params", type=str, default="--retree 2 --maxiterate 0",
                        help="Parameters for mafft, wrap in a single quote (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    return args

def parse_regions(argument):
    """
    Parse the --region rgument
    returns list of regions
    """
    ret = []
    if not argument:
        return ret
    if os.path.exists(argument):
        with open(argument, 'r') as fh:
            for i in fh:
                chrom, start, end = i.strip().split('\t')[:3]
                start = int(start)
                end = int(end)
                ret.append((chrom, start, end))
    else:
        for i in argument.split(','):
            chrom, start, end = re.split(':|-', i)
            start = int(start)
            end = int(end)
            ret.append((chrom, start, end))
    return ret

def get_reference(fn, output, chrom, start, end):
    """
    Pull a subset of the reference and put into a file
    """
    fasta = pysam.FastaFile(fn)
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
    chrom, start, end = region
    if samples is not None:
        samples = "-s " + ",".join(samples)
    else:
        samples = ""
    cmd = f"""bcftools view -c 1 {samples} -r {chrom}:{start}-{end} {vcf} \
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
    chrom, start, end = region
    cmd = """samtools faidx {ref} {chrom}:{start}-{end} | \
bcftools consensus -H{hap} --sample {sample} --prefix {sample}_{hap}_ {vcf} >> {output}
"""
    for samp in samples:
        for hap in [1, 2]:
            m_cmd = cmd.format(ref=ref,
                               chrom=chrom,
                               start=start,
                               end=end,
                               hap=hap,
                               sample=samp,
                               vcf=vcf,
                               output=output)
            ret = truvari.cmd_exe(m_cmd, pipefail=True)
            if ret.ret_code != 0:
                logging.error("Unable to make haplotype %d for sample %s", hap, samp)
                logging.error(ret.stderr)
                sys.exit(1)

def run_mafft(seq_fn, output, params="--retree 2 --maxiterate 0"):
    """
    Run mafft
    """
    cmd = f"mafft {params} {seq_fn} > {output}"
    ret = truvari.cmd_exe(cmd)
    if ret.ret_code != 0: # this doesn't (ever?) exit non-zero
        logging.error("Unable to run MAFFT on %s", seq_fn)
        logging.error(ret.stderr)
        sys.exit(1)

def phab(base_vcf, reference, output_dir, var_region, buffer=100,
        comp_vcf=None, bSamples=None, cSamples=None,
        mafft_params="--retree 2 --maxiterate 0"):
    """
    Harmonize variants with MSA.

    :param `base_vcf`: VCF file name
    :type `base_vcf`: :class:`str`
    :param `reference`: Reference file name
    :type `reference`: :class:`str`
    :param `output_dir`: Destination for results
    :type `output_dir`: :class:`str`
    :param `var_region`: Tuple of region's (chrom, start, end)
    :type `var_region`: :class:`tuple`
    :param `buffer`: Flanking reference bases to add to haplotypes
    :type `buffer`: :class:`int`
    :param `comp_vcf`: VCF file name
    :type `comp_vcf`: :class:`str`
    :param `bSamples`: Samples from `base_vcf` to create haplotypes
    :type `bSamples`: :class:`list`
    :param `cSamples`: Samples from `comp_vcf` to create haplotypes
    :type `cSamples`: :class:`list`
    :param `mafft_params`: Parameters for mafft
    :type `mafft_params`: :class:`str`

    Raises IOError if `output_dir` does not exist
    """
    if not os.path.exists(output_dir):
        raise IOError(f"Directory {output_dir} does not exist")

    buff_region = (var_region[0], var_region[1] - buffer, var_region[2] + buffer)

    if bSamples is None:
        bSamples = list(pysam.VariantFile(base_vcf).header.samples)
    if comp_vcf and cSamples is None:
        cSamples = list(pysam.VariantFile(comp_vcf).header.samples)

    sequences = os.path.join(output_dir, "haps.fa")
    get_reference(reference, sequences, *buff_region)

    subset_base_vcf = os.path.join(output_dir, "base.vcf.gz")
    pull_variants(base_vcf, var_region, subset_base_vcf, reference, bSamples)
    build_consensus(subset_base_vcf, reference, buff_region, sequences, bSamples)

    # if add-to, we put this elsewhere and revisit for now I'll disable that feature
    # because it assumes we already have an MSA
    if comp_vcf is not None:
        subset_comp_vcf = os.path.join(output_dir, "comp.vcf.gz")
        pull_variants(comp_vcf, var_region, subset_comp_vcf, reference, cSamples)
        build_consensus(subset_comp_vcf, reference, buff_region, sequences, cSamples)

    msa_output = os.path.join(output_dir, "msa.fa")
    run_mafft(sequences, msa_output, mafft_params)

    output_vcf = os.path.join(output_dir, "output.vcf")
    vcf = pysam.VariantFile(base_vcf)
    n_header = '##fileformat=VCFv4.1\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    n_header += str(vcf.header.contigs[var_region[0]].header_record)

    with open(output_vcf, 'w') as fout:
        fout.write(truvari.msa2vcf(msa_output, n_header))

    truvari.compress_index_vcf(output_vcf)

def check_requirements():
    """
    ensure external programs are in PATH
    """
    check_fail = False
    for prog in ["bcftools", "vcf-sort", "bgzip", "tabix", "samtools", "mafft"]:
        if not shutil.which(prog):
            logging.error("Unable to find `%s` in PATH", prog)
            check_fail = True
    return check_fail

def phab_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)
    truvari.setup_logging(args.debug)
    if check_requirements():
        sys.exit(1)
    #if check_params(args):
        # mainly just if files exist
        #pass

    if args.bSamples is None:
        args.bSamples = list(pysam.VariantFile(args.base).header.samples)
    else:
        args.bSamples = args.bSamples.split(',')

    if args.comp:
        if args.cSamples is None:
            args.cSamples = list(pysam.VariantFile(args.comp).header.samples)
        else:
            args.cSamples = args.cSamples.split(',')



    # Now I can multiprocess this
    all_regions = parse_regions(args.region)
    for region in all_regions:
        # Make sub-directory if there are multiple regions
        m_output = args.output if len(all_regions) == 1 \
                               else os.path.join(args.output, f"{region[0]}:{region[1]}-{region[2]}")
        try:
            os.makedirs(m_output)
        except FileExistsError:
            logging.error("Output directory %s exists.", args.output)
            sys.exit(1)
        phab(args.base, args.reference, m_output, region, args.buffer,
            args.comp, args.bSamples, args.cSamples, args.mafft_params)

    logging.info("Finished phab")
