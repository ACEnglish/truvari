"""
Wrapper around MAFFT and other tools to harmonize phased variants
"""
import os
import re
import sys
import shutil
import logging
import argparse
import multiprocessing
from functools import partial
from io import BytesIO, StringIO
from collections import defaultdict

import pysam
import pyabpoa
from pysam import samtools
from intervaltree import IntervalTree
from pywfa.align import WavefrontAligner
import truvari

DEFAULT_MAFFT_PARAM = "--auto --thread 1"


def parse_regions(argument):
    """
    Parse the --region argument
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
            except Exception:  # pylint: disable=broad-except
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
    n_reg = 0
    with open(out_file_name, 'w') as fout:
        for chrom in sorted(m_dict.keys()):
            intvs = IntervalTree.from_tuples(m_dict[chrom])
            intvs.merge_overlaps()
            for i in sorted(intvs):
                fout.write(f"{chrom}:{i.begin}-{i.end}\n")
                n_reg += 1
    if n_reg == 0:
        logging.critical("No regions to be refined. Exiting")
        sys.exit(0)
    return out_file_name


def extract_reference(reg_fn, ref_fn):
    """
    Pull the reference
    """
    out_fn = truvari.make_temp_filename(suffix='.fa')
    with open(out_fn, 'w') as fout:
        fout.write(samtools.faidx(ref_fn, "-r", reg_fn))
    return out_fn


def incorporate(consensus_sequence, entry, correction):
    """
    Incorporate a variant into a haplotype returning the new correction field
    """
    if entry.alts[0] == '*':
        return correction

    ref_len = len(entry.ref)
    alt_len = len(entry.alts[0]) if entry.alts else 0
    # Need to check it doesn't overlap previous position
    position = entry.pos + correction
    consensus_sequence[position:position + ref_len] = list(entry.alts[0])
    return correction + (alt_len - ref_len)


def make_haplotypes(sequence, entries, o_samp, ref, start, sample):
    """
    Given a reference sequence, set of entries to incorporate, sample name, reference key, and reference start position
    Make the two haplotypes
    """
    haps = (list(sequence), list(sequence))
    correction = [-start, -start]
    for entry in entries:
        if entry.samples[sample]['GT'][0] == 1:
            correction[0] = incorporate(haps[0], entry, correction[0])
        if len(entry.samples[sample]['GT']) > 1 and entry.samples[sample]['GT'][1] == 1:
            correction[1] = incorporate(haps[1], entry, correction[1])
    return (f">{o_samp}_1_{ref}\n{''.join(haps[0])}\n"
            f">{o_samp}_2_{ref}\n{''.join(haps[1])}\n").encode()


#pylint: disable=too-many-locals
def make_consensus(data, ref_fn, passonly=True, max_size=50000):
    """
    Creates consensus sequence from variants
    """
    vcf_fn, sample, prefix = data
    reference = pysam.FastaFile(ref_fn)
    vcf = truvari.VariantFile(vcf_fn)
    o_samp = 'p:' + sample if prefix else sample
    ret = {}

    tree = defaultdict(IntervalTree)
    for ref in list(reference.references):
        chrom, start, end = re.split(':|-', ref)
        start = int(start)
        end = int(end)
        tree[chrom].addi(start, end+1, data=ref)

    vcf_i = vcf.fetch_regions(tree, with_region=True)

    # Don't make haplotypes of non-sequence resolved, non-pass (sometimes), too big variants
    # Note: max_size == -1 now that we're not capping by variant size.
    # pylint: disable=unnecessary-lambda-assignment
    entry_filter = lambda e: \
           e.is_resolved() \
           and (not passonly or not e.is_filtered()) \
           and (max_size == -1 or e.var_size() <= max_size)
    # pylint: enable=unnecessary-lambda-assignment

    cur_key = None
    cur_entries = []
    for entry, key in vcf_i:
        if cur_key is None:
            cur_key = key

        if key != cur_key:
            ref = f"{cur_key[0]}:{cur_key[1].begin}-{cur_key[1].end - 1}"
            ref_seq = reference.fetch(ref)
            ret[ref] = make_haplotypes(ref_seq, filter(entry_filter, cur_entries), o_samp,
                                       ref, cur_key[1].begin, sample)
            cur_key = key
            cur_entries = []
        cur_entries.append(entry)

    if cur_entries:
        ref = f"{cur_key[0]}:{cur_key[1].begin}-{cur_key[1].end - 1}"
        ref_seq = reference.fetch(ref)
        ret[ref] = make_haplotypes(ref_seq, cur_entries, o_samp,
                                   ref, cur_key[1].begin, sample)

    return ret
#pylint: enable=too-many-locals


def make_haplotype_jobs(base_vcf, bSamples=None, comp_vcf=None, cSamples=None, prefix_comp=False):
    """
    Sets up sample parameters for extract haplotypes
    returns list of tuple of (VCF, sample_name, should_prefix) and the sample names
    to expect in the haplotypes' fasta
    """
    ret = []
    if bSamples is None:
        bSamples = list(truvari.VariantFile(base_vcf).header.samples)
    ret.extend([(base_vcf, samp, False) for samp in bSamples])

    if comp_vcf:
        if cSamples is None:
            cSamples = list(truvari.VariantFile(comp_vcf).header.samples)
        ret.extend([(comp_vcf, samp, prefix_comp or samp in bSamples)
                    for samp in cSamples])

    samp_names = sorted({('p:' if prefix else '') + samp
                         for _, samp, prefix in ret})
    return ret, samp_names


def fasta_reader(fa_str):
    """
    Parses a fasta file as a string and yields tuples of (location, entry)
    """
    cur_name = None
    cur_entry = StringIO()
    for i in fa_str.split('\n'):
        if not i.startswith(">"):
            cur_entry.write(i)
            continue
        if cur_name is not None:
            cur_entry.seek(0)
            yield cur_name, cur_entry.read()
        cur_name = i[1:]
        cur_entry = StringIO()
    cur_entry.seek(0)
    yield cur_name, cur_entry.read()


def collect_haplotypes(ref_haps_fn, hap_jobs, threads, passonly=True, max_size=50000):
    """
    Calls extract haplotypes for every hap_job and the reference sequence
    """
    all_haps = defaultdict(BytesIO)
    m_ref = pysam.FastaFile(ref_haps_fn)
    with multiprocessing.Pool(threads, maxtasksperchild=1000) as pool:
        to_call = partial(make_consensus, ref_fn=ref_haps_fn, passonly=passonly, max_size=max_size)
        # Keep imap because determinism
        for pos, haplotypes in enumerate(pool.imap(to_call, hap_jobs)):
            for location, fasta_entry in haplotypes.items():
                cur = all_haps[location]
                cur.write(fasta_entry)
                if pos == len(hap_jobs) - 1: # Write reference for this location at the end
                    cur.write(f">ref_{location}\n{m_ref[location]}\n".encode())
                    cur.seek(0)
                    all_haps[location] = cur.read()
        pool.close()
        pool.join()
    return all_haps.values()


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
    return "".join(ref), "".join(seq)


def run_wfa(seq_bytes):
    """
    Align haplotypes independently with WFA
    Much faster than mafft, but may be less accurate at finding parsimonous representations
    """
    fasta = dict(fasta_reader(seq_bytes.decode()))
    ref_key = [_ for _ in fasta.keys() if _.startswith("ref_")][0]
    reference = fasta[ref_key]
    aligner = WavefrontAligner(reference, span="end-to-end",
                               heuristic="adaptive")
    for haplotype in fasta:
        if haplotype == ref_key:
            continue
        seq = fasta[haplotype]
        aligner.wavefront_align(seq)
        fasta[haplotype] = expand_cigar(seq, reference, aligner.cigartuples)
    return truvari.msa2vcf(fasta)


def run_mafft(seq_bytes, params=DEFAULT_MAFFT_PARAM):
    """
    Run mafft on the provided sequences provided as a bytestr and return msa2vcf lines
    """
    # Dev only method for overwriting test answers - don't use until code is good
    dev_name = None
    if "PHAB_WRITE_MAFFT" in os.environ and os.environ["PHAB_WRITE_MAFFT"] == "1":
        import hashlib  # pylint: disable=import-outside-toplevel
        dev_name = hashlib.md5(seq_bytes, usedforsecurity=False).hexdigest()

    ret = truvari.cmd_exe(f"mafft {params} -", stdin=seq_bytes)
    if ret.ret_code != 0:
        logging.error("Unable to run MAFFT")
        logging.error("stderr: %s", ret.stderr)
        return ""

    if dev_name is not None:
        with open("repo_utils/test_files/external/fake_mafft/lookup/fm_" + dev_name + ".msa", 'w') as fout:
            fout.write(ret.stdout)

    fasta = dict(fasta_reader(ret.stdout))
    return truvari.msa2vcf(fasta)


def run_poa(seq_bytes):
    """
    Run partial order alignment to create msa
    """
    parts = []
    for k, v in fasta_reader(seq_bytes.decode()):
        parts.append((len(v), v, k))
    parts.sort(reverse=True)
    _, seqs, names = zip(*parts)
    aligner = pyabpoa.msa_aligner()
    aln_result = aligner.msa(seqs, False, True)
    return truvari.msa2vcf(dict(zip(names, aln_result.msa_seq)))


# pylint: disable=too-many-arguments, too-many-locals
# This is just how many arguments it takes
def phab(var_regions, base_vcf, ref_fn, output_fn, bSamples=None, buffer=100,
         mafft_params=DEFAULT_MAFFT_PARAM, comp_vcf=None, cSamples=None,
         prefix_comp=False, threads=1, method="mafft", passonly=True, max_size=50000):
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
    if max_size == -1 or max_size >= 50000:
        logging.warning("Harmonizing variants ≥50kbp is not recommended")

    logging.info("Preparing regions")
    region_fn = merged_region_file(var_regions, buffer)

    logging.info("Extracting haplotypes")
    ref_haps_fn = extract_reference(region_fn, ref_fn)
    hap_jobs, samp_names = make_haplotype_jobs(base_vcf, bSamples,
                                               comp_vcf, cSamples,
                                               prefix_comp)
    haplotypes = collect_haplotypes(ref_haps_fn, hap_jobs, threads, passonly, max_size)

    logging.info("Harmonizing variants")
    if method == "mafft":
        align_method = partial(run_mafft, params=mafft_params)
    elif method == "wfa":
        align_method = run_wfa
    else:
        align_method = run_poa

    with open(output_fn[:-len(".gz")], 'w') as fout:
        fout.write(('##fileformat=VCFv4.1\n'
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'))
        for ctg in truvari.VariantFile(base_vcf).header.contigs.values():
            fout.write(str(ctg.header_record))
        fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        fout.write("\t".join(samp_names) + "\n")

        with multiprocessing.Pool(threads, maxtasksperchild=1000) as pool:
            for result in pool.imap_unordered(align_method, haplotypes):
                fout.write(result)
            pool.close()
            pool.join()

    truvari.compress_index_vcf(output_fn[:-len(".gz")], output_fn)
# pylint: enable=too-many-arguments, too-many-locals

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
    parser.add_argument("-f", "--reference", type=str, required=True,
                        help="Reference")
    parser.add_argument("-o", "--output", type=str, default="phab_out.vcf.gz",
                        help="Output VCF")
    parser.add_argument("--buffer", type=int, default=100,
                        help="Number of reference bases before/after region to add to MSA (%(default)s)")
    parser.add_argument("--align", type=str, choices=["mafft", "wfa", "poa"], default="poa",
                        help="Alignment method (%(default)s)")
    parser.add_argument("-m", "--mafft-params", type=str, default=DEFAULT_MAFFT_PARAM,
                        help="Parameters for mafft, wrap in a single quote (%(default)s)")
    parser.add_argument("--bSamples", type=str, default=None,
                        help="Subset of samples to MSA from base-VCF")
    parser.add_argument("-c", "--comp", type=str, default=None,
                        help="Comparison vcf to MSA")
    parser.add_argument("--cSamples", type=str, default=None,
                        help="Subset of samples to MSA from comp-VCF")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads (%(default)s)")
    parser.add_argument("--maxsize", type=int, default=50000,
                        help="Maximum size of variant to incorporate into haplotypes (%(default)s)")
    parser.add_argument("--passonly", action="store_true",
                        help="Only incorporate passing variants (%(default)s)")
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
        if not shutil.which("mafft"):
            logging.error("Unable to find mafft in PATH")
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
    if os.path.exists(args.output):
        logging.error("Output file already exists")
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
    if args.comp is not None and not truvari.check_vcf_index(args.comp):
        logging.error(
            "Comparison vcf %s must be indexed", args.comp)
        check_fail = True
    if not args.base.endswith(".gz"):
        logging.error(
            "Base vcf %s does not end with .gz. Must be bgzip'd", args.base)
        check_fail = True
    if not truvari.check_vcf_index(args.base):
        logging.error(
            "Base vcf %s must be indexed", args.base)
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
         args.mafft_params, args.comp, args.cSamples, threads=args.threads, method=args.align,
         passonly=args.passonly, max_size=args.maxsize)

    logging.info("Finished phab")
