"""
Wrapper around MAFFT and other tools to perform an MSA comparison of phased variants
"""
import os
import re
import sys
import shutil
import logging
import argparse
import multiprocessing
from io import BytesIO
from functools import partial
from collections import defaultdict

import pysam
from pysam import bcftools, samtools
from intervaltree import IntervalTree
import truvari

try:
    from pywfa.align import WavefrontAligner
    wfa_available=True
except ImportError:
    # Can't load 
    wfa_available=False

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
            try:
                chrom, start, end = i.strip().split('\t')[:3]
                start = int(start)
                end = int(end)
            except Exception: #pylint: disable=broad-except
                logging.error("Unable to parse bed line %s", i)
                sys.exit(1)
            ret.append((chrom, start, end))
    else:
        for i in argument.split(','):
            try:
                chrom, start, end = re.split(':|-', i)
                start = int(start)
                end = int(end)
            except ValueError:
                logging.error("Unable to parse region line %s", i)
                sys.exit(1)
            ret.append((chrom, start, end))
    return ret

def merged_region_file(regions, buff=100):
    """
    Write a file for extracting reference regions with samtools
    returns the name of the temporary file made
    """
    m_dict = defaultdict(list)
    for i in regions:
        m_dict[i[0]].append((max(0, i[1] - buff), i[2] + buff))

    out_file_name = truvari.make_temp_filename()
    with open(out_file_name, 'w') as fout:
        for chrom in m_dict:
            intvs = IntervalTree.from_tuples(m_dict[chrom])
            intvs.merge_overlaps()
            for i in sorted(intvs):
                fout.write(f"{chrom}:{i.begin}-{i.end}\n")
    return out_file_name

def extract_reference(reg_fn, ref_fn):
    """
    Pull the reference
    """
    out_fn = truvari.make_temp_filename(suffix='.fa')
    with open(out_fn, 'w') as fout:
        fout.write(samtools.faidx(ref_fn, "-r", reg_fn))
    return out_fn

def make_haplotype_jobs(base_vcf, bSamples, comp_vcf, cSamples, prefix_comp):
    """
    Sets up sample parameters for extract haplotypes
    returns list of tuple of (VCF, sample_name, should_prefix) and the sample names
    to expect in the haplotypes' fasta
    """
    ret = []
    if bSamples is None:
        ret.extend([(base_vcf, samp, False)
                    for samp in list(pysam.VariantFile(base_vcf).header.samples)])
    if comp_vcf and cSamples is None:
        bsamps = list(pysam.VariantFile(base_vcf).header.samples)
        ret.extend([(comp_vcf, samp, prefix_comp or samp in bsamps)
                    for samp in list(pysam.VariantFile(comp_vcf).header.samples)])
    samp_names = [('p:' if prefix else '') + samp for _, samp, prefix in ret]
    samp_names.sort()
    return ret, samp_names

def fasta_reader(fa_str, name_entries=True):
    """
    Parses a fasta file as a string and yields tuples of (location, entry)
    if name_entries, the entry names are written to the value and location is honored
    """
    cur_name = None
    cur_entry = BytesIO()
    for i in fa_str.split('\n'):
        if not i.startswith(">"):
            cur_entry.write(i.encode())
            continue
        if cur_name is not None:
            cur_entry.write(b'\n')
            cur_entry.seek(0)
            yield cur_name, cur_entry.read()
        cur_name = i[1:]
        cur_entry = BytesIO()
        if name_entries:
            cur_name = cur_name.split('_')[-1]
            cur_entry.write((i + '\n').encode())

    cur_entry.write(b'\n')
    cur_entry.seek(0)
    yield cur_name, cur_entry.read()

def extract_haplotypes(data, ref_fn):
    """
    Given a data tuple of VCF, sample name, and prefix bool
    Call bcftools consensus
    Returns list of tuples [location, fastaentry] for every haplotype, so locations are duplicated
    """
    vcf_fn, sample, prefix = data
    prefix = 'p:' if prefix else ''
    ret = []
    cmd = "-H{{hap}} --sample {sample} --prefix {prefix}{sample}_{{hap}}_ -f {ref} {vcf}"
    cmd = cmd.format(sample=sample, prefix = prefix, ref=ref_fn, vcf=vcf_fn)
    for i in [1, 2]:
        ret.extend(fasta_reader(bcftools.consensus(*cmd.format(hap=i).split(' '))))
    return ret

def collect_haplotypes(ref_haps_fn, hap_jobs, threads):
    """
    Calls extract haplotypes for every hap_job
    """
    all_haps = defaultdict(BytesIO)
    with multiprocessing.Pool(threads, maxtasksperchild=1) as pool:
        for haplotype in pool.imap_unordered(partial(extract_haplotypes, ref_fn=ref_haps_fn),
                                             hap_jobs):
            for location, fasta_entry in haplotype:
                all_haps[location].write(fasta_entry)

    m_ref = pysam.FastaFile(ref_haps_fn)
    for location in list(m_ref.references):
        all_haps[location].write(f">ref_{location}\n".encode())
        all_haps[location].write(m_ref[location].encode())
        all_haps[location].write(b'\n')
        all_haps[location].seek(0)
        yield all_haps[location].read()

def expand_cigar(seq, ref, cigar):
    """
    Put gaps '-' in sequence where needed
    """
    seq_pos = 0
    seq = list(seq)
    ref_pos = 0
    ref = list(ref)
    for code, span in cigar:
        if code == 2:
            seq.insert(seq_pos, '-' * span)
            seq_pos += 1
            ref_pos += span
        elif code == 1:
            ref.insert(ref_pos, '-' * span)
            seq_pos += span
            ref_pos += 1
        else:
            seq_pos += span
            ref_pos += span
    print("".join(seq))
    return "".join(seq), "".join(ref)

def wfa_to_vars(seq_bytes):
    """
    """
    fasta = {k:v.decode() for k,v in fasta_reader(seq_bytes.decode(), name_entries=False)}
    ref_key = [_ for _ in fasta.keys() if _.startswith("ref_")][0]
    reference = fasta[ref_key]
    # Get the ref sequence
    aligner = WavefrontAligner(scope="score",
                               distance="affine2p",
                               span="end-to-end",
                               heuristic="adaptive")
    for haplotype in fasta:
        if haplotype == ref_key:
            continue
        seq = fasta[haplotype]
        aligner.wavefront_align(reference, seq)
        assert aligner.status == 0
        fasta[haplotype] = expand_cigar(seq, reference, aligner.cigartuples)
    return truvari.msa2vcf(fasta)

def edlib_to_vars(seq_bytes):
    import edlib
    lookup = defaultdict(int)
    lookup.update({"I":1, "D": 2})
    def to_tuple(cig):
        cur_pos = 0
        i_pos = 0
        while i_pos < len(cig):
            if cig[i_pos].isdigit():
                i_pos += 1
            else:
                s = (lookup[cig[i_pos]], int(cig[cur_pos:i_pos]))
                yield s
                i_pos += 1
                cur_pos = i_pos
        
        
    fasta = {k:v.decode() for k,v in fasta_reader(seq_bytes.decode(), name_entries=False)}
    ref_key = [_ for _ in fasta.keys() if _.startswith("ref_")][0]
    reference = fasta[ref_key]
    # Get the ref sequence

    for haplotype in fasta:
        if haplotype == ref_key:
            continue
        seq = fasta[haplotype]
        result = edlib.align(reference, seq, mode = "HW", task = "path")
        fasta[haplotype] = expand_cigar(seq, reference, to_tuple(result["cigar"]))
    return truvari.msa2vcf(fasta)

    
DEFAULT_MAFFT_PARAM="--auto --thread 1"
def mafft_to_vars(seq_bytes, params=DEFAULT_MAFFT_PARAM):
    """
    Run mafft on the provided sequences provided as a bytestr and return msa2vcf lines
    """
    # Dev only method for overwriting test answers - don't use until code is good
    dev_name = None
    if "PHAB_WRITE_MAFFT" in os.environ and os.environ["PHAB_WRITE_MAFFT"] == "1":
        import hashlib # pylint: disable=import-outside-toplevel
        dev_name = hashlib.md5(seq_bytes).hexdigest()

    cmd = f"mafft --quiet {params} -"
    ret = truvari.cmd_exe(cmd, stdin=seq_bytes)
    if ret.ret_code != 0:
        logging.error("Unable to run MAFFT")
        logging.error(ret.stderr)
        return ""

    if dev_name is not None:
        with open("repo_utils/test_files/external/fake_mafft/lookup/fm_" + dev_name + ".msa", 'w') as fout:
            fout.write(ret.stdout)

    fasta = {}
    for name, sequence in fasta_reader(ret.stdout, name_entries=False):
        fasta[name] = sequence.decode()
    return truvari.msa2vcf(fasta)

def harmonize_variants(harm_jobs, mafft_params, base_vcf, samp_names, output_fn, threads, method="mafft"):
    """
    Parallel processing of variants to harmonize. Writes to output
    """
    # output_fn, base_vcf, samp_names, threads, harm_jobs
    align_method = partial(mafft_to_vars, params=mafft_params)
    if method == "wfa":
        align_method = wfa_to_vars
    elif method == "edlib":
        align_method = edlib_to_vars
    with open(output_fn[:-len(".gz")], 'w') as fout:
        # Write header
        fout.write('##fileformat=VCFv4.1\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        for ctg in pysam.VariantFile(base_vcf).header.contigs.values():
            fout.write(str(ctg.header_record))
        fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        fout.write("\t".join(samp_names) + "\n")
        # Write mafft
        with multiprocessing.Pool(threads) as pool:
            for result in pool.imap_unordered(align_method, harm_jobs):
                fout.write(result)
    truvari.compress_index_vcf(output_fn[:-len(".gz")], output_fn)


#pylint: disable=too-many-arguments
# This is just how many arguments it takes
def phab(var_regions, base_vcf, ref_fn, output_fn, bSamples=None, buffer=100,
         mafft_params=DEFAULT_MAFFT_PARAM, comp_vcf=None, cSamples=None,
         prefix_comp=False, threads=1, method="mafft"):
    """
    Harmonize variants with MSA. Runs on a set of regions given files to create
    haplotypes and writes results to a compressed/indexed VCF

    :param `var_regions`: List of tuples of region's (chrom, start, end)
    :type `var_regions`: :class:`list`
    :param `base_vcf`: VCF file name
    :type `base_vcf`: :class:`str`
    :param `ref_fn`: Reference file name
    :type `ref_fn`: :class:`str`
    :param `output_fn`: Output VCF.gz
    :param `bSamples`: Samples from `base_vcf` to create haplotypes
    :type `bSamples`: :class:`list`
    :param `buffer`: Flanking reference bases to added to regions
    :type `buffer`: :class:`int`
    :param `mafft_params`: Parameters for mafft
    :type `mafft_params`: :class:`str`
    :param `comp_vcf`: VCF file name
    :type `comp_vcf`: :class:`str`
    :param `cSamples`: Samples from `comp_vcf` to create haplotypes
    :type `cSamples`: :class:`list`
    :param `prefix_comp`: Ensure unique sample names by prefixing comp samples
    :type `prefix_comp`: :class:`bool`
    :param `threads`: Number of threads to use
    :type `threads`: :class:`int`
    :param `method`: Alignment method to use mafft or wfa
    :type `method`: :class:`str`
    """
    logging.info("Preparing regions")
    region_fn = merged_region_file(var_regions, buffer)

    logging.info("Extracting haplotypes")
    ref_haps_fn = extract_reference(region_fn, ref_fn)
    hap_jobs, samp_names = make_haplotype_jobs(base_vcf, bSamples, comp_vcf, cSamples, prefix_comp)
    harm_jobs = collect_haplotypes(ref_haps_fn, hap_jobs, threads)

    logging.info("Harmonizing variants")
    harmonize_variants(harm_jobs, mafft_params, base_vcf, samp_names, output_fn, threads, method)
#pylint: enable=too-many-arguments

######
# UI #
######
def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="phab", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-r", "--region", type=str, required=True,
                        help="Bed filename or comma-separated list of chrom:start-end regions to process")
    parser.add_argument("-b", "--base", type=str, required=True,
                        help="Baseline vcf to MSA")
    parser.add_argument("-c", "--comp", type=str, default=None,
                        help="Comparison vcf to MSA")
    parser.add_argument("-f", "--reference", type=str, required=True,
                        help="Reference")
    parser.add_argument("--buffer", type=int, default=100,
                        help="Number of reference bases before/after region to add to MSA (%(default)s)")
    parser.add_argument("-o", "--output", type=str, default="phab_out.vcf.gz",
                        help="Output VCF")
    parser.add_argument("--bSamples", type=str, default=None,
                        help="Subset of samples to MSA from base-VCF")
    parser.add_argument("--cSamples", type=str, default=None,
                        help="Subset of samples to MSA from comp-VCF")
    parser.add_argument("-m", "--mafft-params", type=str, default=DEFAULT_MAFFT_PARAM,
                        help="Parameters for mafft, wrap in a single quote (%(default)s)")
    parser.add_argument("--align", type=str, choices=["mafft", "edlib", "wfa"], default="mafft",
                        help="Alignment method accurate (mafft) or fast (wfa) (%(default)s)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    return args

def check_requirements(align):
    """
    ensure external programs are in PATH
    """
    check_fail = False
    if align == "mafft":
        if not shutil.which(prog):
            logging.error("Unable to find `%s` in PATH", prog)
            check_fail = True
    elif align == "wfa" and not wfa_available:
        logging.error("pywfa could not be imported")
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

def phab_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)
    truvari.setup_logging(args.debug, show_version=True)
    if check_requirements(args.align) or check_params(args):
        logging.error("Couldn't run Truvari. Please fix parameters\n")
        sys.exit(100)

    if args.bSamples is not None:
        args.bSamples = args.bSamples.split(',')

    if args.comp and args.cSamples is not None:
        args.cSamples = args.cSamples.split(',')

    all_regions = parse_regions(args.region)
    phab(all_regions, args.base, args.reference, args.output, args.bSamples, args.buffer,
         args.mafft_params, args.comp, args.cSamples, threads=args.threads, method=args.align)

    logging.info("Finished phab")
