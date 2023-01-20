"""
Wrapper around MAFFT and other tools to perform an MSA comparison of phased variants
"""
import os
import re
import sys
import glob
import shutil
import logging
import argparse
import multiprocessing

import pysam
from pysam import bcftools
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
                        help="Number of reference bases before/after region to add to MSA (%(default)s)")
    #parser.add_argument("--add-to", action='store_true',
    #                    help="Build the baseMSA independentally of the compMSA, then add the comp")
    parser.add_argument("-o", "--output", type=str, default="phab_out.vcf.gz",
                        help="Output VCF")
    parser.add_argument("-k", "--keep-parts", type=str, default=None,
                        help="Directory to save intermediate region results")
    parser.add_argument("--bSamples", type=str, default=None,
                        help="Subset of samples to MSA from base-VCF")
    parser.add_argument("--cSamples", type=str, default=None,
                        help="Subset of samples to MSA from comp-VCF")
    parser.add_argument("-m", "--mafft-params", type=str, default="--auto",
                        help="Parameters for mafft, wrap in a single quote (%(default)s)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads (%(default)s)")
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
        for i in truvari.opt_gz_open(argument):
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
    cnt = 0
    with pysam.VariantFile(output) as fh:
        for _ in fh:
            cnt += 1
    return cnt

def build_consensus(vcf, ref, region, output, samples=None, prefix_name=False):
    """
    Make the consensus sequence - appends to output
    """
    chrom, start, end = region
    prefix = 'p:' if prefix_name else ''
    cmd = """samtools faidx {ref} {chrom}:{start}-{end} | \
bcftools consensus -H{hap} --sample {sample} --prefix {prefix}{sample}_{hap}_ {vcf} >> {output}
"""
    for samp in samples:
        for hap in [1, 2]:
            m_cmd = cmd.format(ref=ref,
                               chrom=chrom,
                               start=start,
                               end=end,
                               hap=hap,
                               sample=samp,
                               prefix=prefix,
                               vcf=vcf,
                               output=output)
            ret = truvari.cmd_exe(m_cmd, pipefail=True)
            if ret.ret_code != 0:
                logging.error("Unable to make haplotype %d for sample %s", hap, samp)
                logging.error(ret.stderr)
                sys.exit(1)

def run_mafft(seq_fn, output, params="--auto"):
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
        mafft_params="--auto", prefix_comp=False):
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
    :param `prefix_comp`: Ensure unique sample names by prefixing comp samples
    :type `prefix_comp`: :class:`bool`

    Raises IOError if `output_dir` does not exist
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    buff_region = (var_region[0], var_region[1] - buffer, var_region[2] + buffer)

    if bSamples is None:
        bSamples = list(pysam.VariantFile(base_vcf).header.samples)
    if comp_vcf and cSamples is None:
        cSamples = list(pysam.VariantFile(comp_vcf).header.samples)

    # Make sure there are variants first
    subset_base_vcf = os.path.join(output_dir, "base.vcf.gz")
    num_vars = pull_variants(base_vcf, var_region, subset_base_vcf, reference, bSamples)
    if comp_vcf is not None:
        subset_comp_vcf = os.path.join(output_dir, "comp.vcf.gz")
        num_vars += pull_variants(comp_vcf, var_region, subset_comp_vcf, reference, cSamples)

    if num_vars == 0:
        return

    sequences = os.path.join(output_dir, "haps.fa")
    get_reference(reference, sequences, *buff_region)

    build_consensus(subset_base_vcf, reference, buff_region, sequences, bSamples)

    # if add-to, we put this elsewhere and revisit for now I'll disable that feature
    # because it assumes we already have an MSA
    if comp_vcf is not None:
        build_consensus(subset_comp_vcf, reference, buff_region, sequences, cSamples, prefix_comp)

    msa_output = os.path.join(output_dir, "msa.fa")
    run_mafft(sequences, msa_output, mafft_params)

    output_vcf = os.path.join(output_dir, "output.vcf")
    vcf = pysam.VariantFile(base_vcf)
    n_header = '##fileformat=VCFv4.1\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    n_header += str(vcf.header.contigs[var_region[0]].header_record)

    with open(output_vcf, 'w') as fout:
        try:
            fout.write(truvari.msa2vcf(msa_output, n_header))
        except RuntimeWarning:
            return

    truvari.compress_index_vcf(output_vcf)

def consolidate_phab_vcfs(phab_dir, out_vcf):
    """
    Consolidate all the phab output VCFs in a directory and write to compressed indexed vcf
    """
    def concat(file_names):
        tmp_name = truvari.make_temp_filename(suffix=".vcf.gz")
        files = truvari.make_temp_filename(suffix=".txt")
        with open(files, 'w') as fout:
            fout.write('\n'.join(file_names))
        bcftools.concat("--no-version", "-O", "z", "-a", "-f", files, "-o", tmp_name, catch_stdout=False)
        pysam.tabix_index(tmp_name, preset='vcf')
        return tmp_name

    in_files = glob.glob(os.path.join(phab_dir, "*", "output.vcf.gz"))
    in_files.sort()
    while len(in_files) > 1:
        tmp_names = []
        # bcftools gets weird with more than 1,010 files
        for i in range(0, len(in_files), 1010):
            tmp_names.append(concat(in_files[i:i+1010]))
        in_files = tmp_names
    shutil.move(in_files[0], out_vcf)
    shutil.move(in_files[0] + '.tbi', out_vcf + '.tbi')

def check_requirements():
    """
    ensure external programs are in PATH
    """
    check_fail = False
    for prog in ["bcftools", "bgzip", "tabix", "samtools", "mafft"]:
        if not shutil.which(prog):
            logging.error("Unable to find `%s` in PATH", prog)
            check_fail = True
    return check_fail

def check_params(args):
    """
    Ensure files are okay to use
    """
    check_fail = False
    if not args.output.endswith(".vcf.gz"):
        logging.error("Output file must be a '.vcf.gz', got %s", args.output)
        check_fail = True
    if args.keep_parts and os.path.isdir(args.keep_parts):
        logging.error("Output directory '%s' already exists", args.keep_parts)
        check_fail = True
    if args.comp is not None and not os.path.exists(args.comp):
        logging.error("File %s does not exist", args.comp)
        check_fail = True
    if not os.path.exists(args.base):
        logging.error("File %s does not exist", args.base)
        check_fail = True
    if args.comp is not None and not args.comp.endswith(".gz"):
        logging.error(
            "Comparison vcf %s does not end with .gz. Must be bgzip'd", args.comp)
        check_fail = True
    if args.comp is not None and not os.path.exists(args.comp + '.tbi'):
        logging.error(
            "Comparison vcf index %s.tbi does not exist. Must be indexed", args.comp)
        check_fail = True
    if not args.base.endswith(".gz"):
        logging.error(
            "Base vcf %s does not end with .gz. Must be bgzip'd", args.base)
        check_fail = True
    if not os.path.exists(args.base + '.tbi'):
        logging.error(
            "Base vcf index %s.tbi does not exist. Must be indexed", args.base)
        check_fail = True
    if not os.path.exists(args.reference):
        logging.error("Reference %s does not exist", args.reference)
        check_fail = True
    return check_fail

def phab_wrapper(job):
    """
    Convience function to call phab (which works no a single region)
    job is a tuple of region, output_dir, and kwargs dict
    """
    try:
        phab(var_region=job[0], output_dir=job[1], **job[2])
    except Exception as e: #pylint: disable=broad-except
        logging.critical("phab failed on %s\n%s", str(job[0]), e)

#pylint: disable=too-many-arguments
# This is just how many arguments it takes
def phab_multi(base_vcf, reference, output_dir, var_regions, buffer=100, comp_vcf=None,
               bSamples=None, cSamples=None, mafft_params="--auto",
               prefix_comp=False, threads=1):
    """
    Run phab on multiple regions
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    params = {"base_vcf": base_vcf,
              "reference": reference,
              "buffer": buffer,
              "comp_vcf": comp_vcf,
              "bSamples": bSamples,
              "cSamples": cSamples,
              "mafft_params": mafft_params,
              "prefix_comp": prefix_comp}
    with multiprocessing.Pool(threads) as pool:
        # build jobs
        jobs = []
        for region in var_regions:
            m_output = os.path.join(output_dir, f"{region[0]}:{region[1]}-{region[2]}")
            jobs.append((region, m_output, params))

        pool.imap_unordered(phab_wrapper, jobs)
        pool.close()
        pool.join()
#pylint: enable=too-many-arguments

def phab_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)
    truvari.setup_logging(args.debug, show_version=True)
    if check_requirements() or check_params(args):
        logging.error("Couldn't run Truvari. Please fix parameters\n")
        sys.exit(100)

    # Setting up the samples is tricky
    if args.bSamples is None:
        args.bSamples = list(pysam.VariantFile(args.base).header.samples)
    else:
        args.bSamples = args.bSamples.split(',')

    if args.comp:
        if args.cSamples is None:
            args.cSamples = list(pysam.VariantFile(args.comp).header.samples)
        else:
            args.cSamples = args.cSamples.split(',')

    prefix_comp = False
    if args.cSamples and set(args.cSamples) & set(args.bSamples):
        logging.warning("--cSamples intersect with --bSamples. Output vcf --comp SAMPLE names will have prefix 'p:'")
        prefix_comp = True

    all_regions = parse_regions(args.region)
    if args.keep_parts is None:
        remove = True # remove at the end
        args.keep_parts = truvari.make_temp_filename()
    else:
        remove = False
    phab_multi(args.base, args.reference, args.keep_parts, all_regions, args.buffer, args.comp,
               args.bSamples, args.cSamples, args.mafft_params, prefix_comp, args.threads)

    consolidate_phab_vcfs(args.keep_parts, args.output)

    if remove:
        logging.info("Cleaning")
        shutil.rmtree(args.keep_parts)

    logging.info("Finished phab")
