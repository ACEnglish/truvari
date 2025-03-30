"""
Wrapper around MAFFT and other tools to harmonize phased variants
"""
import os
import re
import sys
import time
import atexit
import pickle
import shutil
import logging
import argparse
import functools
from io import StringIO
from collections import defaultdict
import multiprocessing
import multiprocessing.shared_memory as shm


import pysam
import pyabpoa
from pysam import samtools
from intervaltree import IntervalTree
from pywfa.align import WavefrontAligner
import truvari

DEFAULT_MAFFT_PARAM = "--auto --thread 1"

##################
# phab utilities #
##################


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


def make_haplotypes(sequence, entries, ref, start, in_sample, out_sample):
    """
    Given a reference sequence, set of entries to incorporate, sample name, reference key, and reference start position
    Make the two haplotypes
    """
    haps = (list(sequence), list(sequence))
    correction = [-start, -start]
    for entry in entries:
        if entry.samples[in_sample]['GT'][0] == 1:
            correction[0] = incorporate(haps[0], entry, correction[0])
        if len(entry.samples[in_sample]['GT']) > 1 and entry.samples[in_sample]['GT'][1] == 1:
            correction[1] = incorporate(haps[1], entry, correction[1])
    return {f"{out_sample}_1_{ref}": ''.join(haps[0]),
            f"{out_sample}_2_{ref}": ''.join(haps[1])}

#####################
# aligner utilities #
#####################


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

def deduplicate_haps(d):
    """
    Deduplicates a dictionary by replacing duplicate values with a single key.
    """
    value_to_key = {}
    dedup_dict = {}
    key_mapping = {}

    for key, value in d.items():
        if value in value_to_key:
            key_mapping[key] = value_to_key[value]
        else:
            dedup_key = key  # Use the first occurrence as the deduplicated key
            value_to_key[value] = dedup_key
            dedup_dict[dedup_key] = value
            key_mapping[key] = dedup_key

    return dedup_dict, key_mapping

def safe_align_method(job, func, dedup=True):
    """
    Wrapper safely calling the alignment method on a PhabJob or dictionary
    Will then return a call to truvari.msa2vcf on the align method's result
    Note, if the MSA fails, this return s a string "ERROR: {e}" instead of raising an exception
    """
    if isinstance(job, PhabJob):
        work = job.build_haplotypes()
    elif isinstance(job, dict):
        work = job
    else:
        raise TypeError(f"Unknown job type {type(job)}")

    if dedup:
        work, key_map = deduplicate_haps(work)

    try:
        msa = func(work)
    except Exception as e:  # pylint: disable=broad-exception-caught
        return f"ERROR: {e}"

    if dedup:
        msa = {key: msa[val] for key,val in key_map.items()}

    return truvari.msa2vcf(msa)

def run_wfa(haplotypes, aligner=None):
    """
    Align haplotypes independently with WFA
    Much faster than mafft, but may be less accurate at finding parsimonous representations
    a pre-build `WavefrontAligner` may be passed
    """
    ref_key = [_ for _ in haplotypes.keys() if _.startswith("ref_")][0]
    reference = haplotypes[ref_key]
    aligner = WavefrontAligner(reference, span="end-to-end",
                                  heuristic="adaptive")
    for haplotype in haplotypes:
        if haplotype == ref_key:
            continue
        seq = haplotypes[haplotype]
        aligner.wavefront_align(seq)
        haplotypes[haplotype] = expand_cigar(
            seq, reference, aligner.cigartuples)
    return haplotypes

def run_mafft(haplotypes, params=DEFAULT_MAFFT_PARAM):
    """
    Run mafft on the provided haplotypes dictionary
    """
    seq_bytes = "".join(
        [f'>{k}\n{v}\n' for k, v in haplotypes.items()]).encode()

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
    return fasta


def run_poa(haplotypes, aligner=None, out_cons=False, out_msa=True):
    """
    Run partial order alignment of haplotypes to create msa
    a pre-build pyabpoa.msa_aligner may be passed
    """
    parts = []
    for k, v in haplotypes.items():
        parts.append((len(v), v, k))
    parts.sort(reverse=True)
    _, seqs, names = zip(*parts)
    if aligner is None:
        aligner = pyabpoa.msa_aligner()
    aln_result = aligner.msa(seqs, out_cons, out_msa)
    return dict(zip(names, aln_result.msa_seq))


class VCFtoHaplotypes():
    """
    Class for holding input VCFs with helpers for fetching/filtering
    variants and writing the output header
    """

    def __init__(self, reference_fn, vcf_fns, samples=None, passonly=True, max_size=50000):
        self.reference_fn = reference_fn
        self.vcf_fns = vcf_fns
        self.passonly = passonly
        self.max_size = max_size
        # Filled in later by `self.set_regions`
        self.ref_haps_fn = None

        self.sample_lookup = {}
        # Counter for unique suffix on output names
        seen_samples = {}
        self.out_samples = []
        for i in vcf_fns:
            self.sample_lookup[i] = []
            for old_sample in list(truvari.VariantFile(i).header.samples):
                # Only grab the requested samples
                if samples is not None and old_sample not in samples:
                    continue
                # --force-sample when needed
                if old_sample in seen_samples:
                    new_sample = old_sample + ':' + str(seen_samples[old_sample])
                    seen_samples[old_sample] += 1
                else:
                    new_sample = old_sample
                    seen_samples[old_sample] = 1
                # For the output vcf
                self.out_samples.append(new_sample)
                # For parsing input vcfs
                self.sample_lookup[i].append((old_sample, new_sample))

        # Must be sorted for the msa2vcf later
        self.out_samples.sort()

    def set_regions(self, regions, buff=100):
        """
        Write a file for extracting reference regions with samtools
        Then extracts those regions.
        Sets the class's x/y variables
        """
        m_dict = defaultdict(list)
        for i in regions:
            m_dict[i[0]].append((max(0, i[1] - buff), i[2] + buff))

        regions_file_name = truvari.make_temp_filename()
        n_reg = 0
        with open(regions_file_name, 'w') as fout:
            for chrom in sorted(m_dict.keys()):
                intvs = IntervalTree.from_tuples(m_dict[chrom])
                intvs.merge_overlaps()
                for i in sorted(intvs):
                    fout.write(f"{chrom}:{i.begin}-{i.end}\n")
                    n_reg += 1
        if n_reg == 0:
            logging.critical("No regions to be refined. Exiting")
            sys.exit(0)

        # Pull sequences
        out_fn = truvari.make_temp_filename(suffix='.fa')
        with open(out_fn, 'w') as fout:
            fout.write(samtools.faidx(
                self.reference_fn, "-r", regions_file_name))
        # Facilitate fetching
        samtools.faidx(out_fn)
        self.ref_haps_fn = out_fn

    def get_haplotypes(self, refname):
        """
        Fetches the variants and builds the haplotypes of all the self.samples_lookup
        returns the fasta dictionary
        """
        chrom, start, end = re.split(':|-', refname)
        start = int(start)
        end = int(end)
        refseq = pysam.FastaFile(self.ref_haps_fn).fetch(refname)

        filt = functools.partial(self.__keep_entry, start=start, end=end)

        ret = {}
        for vcf_fn in self.vcf_fns:
            vcf = truvari.VariantFile(vcf_fn)
            entries = list(filter(filt, vcf.fetch(chrom, start, end)))
            for in_samp, out_samp in self.sample_lookup[vcf_fn]:
                ret.update(make_haplotypes(refseq, entries, refname,
                                           start, in_samp, out_samp))
            vcf.close()

        ret[f"ref_{refname}"] = refseq
        return ret

    def build_all(self):
        """
        Yields each locus' name and haplotypes as a dict
        """
        for name in list(pysam.FastaFile(self.ref_haps_fn).references):
            yield name, self.get_haplotypes(name)

    def __keep_entry(self, e, start, end):
        """
        I feel like this could be handled by truvari v5 api
        """
        return e.is_resolved() \
            and (not self.passonly or not e.is_filtered()) \
            and (self.max_size == -1 or e.var_size() <= self.max_size) \
            and e.within(start, end)

##################
# Infrastructure #
##################

# pylint: disable=too-few-public-methods
class PhabJob():
    """
    This holds all the information needed to run one task through VCFtoHaplotypes
    It is responsible for fetching reference sequence, sending to VCFtoHaplotypes,
    and then passing that consolidated set of sequences to the align_method
    # This could have a .run that is safe_align_method.. but whatever
    """

    def __init__(self, name, mem_vcf_info):
        self.name = name
        self.mem_vcf_info = mem_vcf_info

    def build_haplotypes(self):
        """
        Returns the fasta dict of haplotypes for this job
        This is essentially my collect_haplotypes
        """
        # Attach shared memory
        shm_name, shm_size = self.mem_vcf_info
        existing_shm = shm.SharedMemory(name=shm_name)
        data = bytes(existing_shm.buf[:shm_size])
        vcf_info = pickle.loads(data)

        return vcf_info.get_haplotypes(self.name)
# pylint: enable=too-few-public-methods


def monitored_pool(method, locus_jobs, threads):
    """
    Create a pool of workers, send jobs to a safe_align_method, monitor the results for yielding
    """
    # given timeout in minutes, figure out how many retries are given
    if "PHAB_TIMEOUT" in os.environ:
        MAXRETRIES = int(os.environ["PHAB_TIMEOUT"]) * 12
    else:
        MAXRETRIES = 60
    WAITINTERVAL = 5
    FAILTIME = time.strftime(
        "%H:%M:%S", time.gmtime(MAXRETRIES * WAITINTERVAL))
    n_completed = 0
    prev_completed = 0.05
    n_failed = 0
    with multiprocessing.Pool(threads, maxtasksperchild=1000) as pool:
        results = [pool.apply_async(method, (job,)) for job in locus_jobs]
        pending = {idx: (result, 0)
                   for idx, result in enumerate(results)}  # Tracks retries

        while pending:
            for idx in list(pending.keys()):
                result, retries = pending[idx]

                try:
                    if result.ready():
                        output = result.get()
                        del pending[idx]

                        if isinstance(output, str) and output.startswith("ERROR:"):
                            job = locus_jobs[idx]
                            logging.error(f"{job.name} {output}")
                            n_failed += 1
                        else:
                            n_completed += 1
                            yield output
                    # only penalize jobs that had a reasonable chance to start
                    elif idx < (n_completed + n_failed + threads):
                        if retries < MAXRETRIES:
                            pending[idx] = (result, retries + 1)
                        else:
                            del pending[idx]
                            n_failed += 1
                            job = locus_jobs[idx]
                            logging.error(
                                f"{job.name} ERROR: Timeout after {FAILTIME}")
                except Exception as e:  # pylint: disable=broad-exception-caught
                    # Remove failed job
                    del pending[idx]
                    n_failed += 1
                    job = locus_jobs[idx]
                    logging.error(f"{job.name} ERROR: {e}")

            if pending:
                pct_completed = (n_completed + n_failed) / len(locus_jobs)
                if pct_completed >= prev_completed:
                    logging.info("Completed %d / %d (%d%%) Loci; %d Failed",
                                 n_completed, len(locus_jobs),
                                 pct_completed * 100, n_failed)
                    prev_completed = min(1, pct_completed + 0.05)
                time.sleep(WAITINTERVAL)  # Wait before retrying

    # I should perhaps exit non-zero if too many jobs fail...


def cleanup_shared_memory(shared_info):
    """
    cleans shm
    """
    try:
        shared_info.close()
        shared_info.unlink()
    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Error cleaning up shared memory: {e}")


def run_phab(vcf_info, regions, output_fn, buffer=100,
             align_method=run_poa, dedup=True, threads=1):
    """
    Harmonize variants with MSA. Runs on a set of regions given files to create
    haplotypes and writes results to a compressed/indexed VCF

    :param `vcf_info` : VCF to Haplotype maker
    :type `vcf_info` : :class:`VCFtoHaplotypes`
    :param `regions`: List of tuples of region's (chrom, start, end)
    :type `regions`: :class:`list`
    :param `reference_fn`: Reference file name
    :type `reference_fn`: :class:`str`
    :param `output_fn`: Output VCF.gz
    :type `output_fn`: :class:`str`
    :param `buffer`: Flanking reference bases to added to regions
    :type `buffer`: :class:`int`
    :param `align_method`: Alignment method's function. See `get_align_method`
    :type `method`: function
    :param `threads`: Number of threads to use
    :type `threads`: :class:`int`
    """
    logging.info("Preparing reference/regions")
    vcf_info.set_regions(regions, buffer)

    logging.info("Harmonizing variants")
    # Shared memory for vcf_info to reduce copies/memory
    data = pickle.dumps(vcf_info)
    shared_info = shm.SharedMemory(create=True, size=len(data))
    shared_info.buf[:len(data)] = data
    mem_vcf_info = (shared_info.name, shared_info.size)
    # Cleanup shared memory on exit
    atexit.register(cleanup_shared_memory, shared_info)

    # Now we build all the jobs
    jobs = [PhabJob(name, mem_vcf_info) for name in
            list(pysam.FastaFile(vcf_info.ref_haps_fn).references)]

    to_call = functools.partial(safe_align_method, func=align_method, dedup=dedup)
    with open(output_fn[:-len(".gz")], 'w') as fout:
        fout.write(('##fileformat=VCFv4.1\n'
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'))
        for ctg in truvari.VariantFile(vcf_info.vcf_fns[0]).header.contigs.values():
            fout.write(str(ctg.header_record))
        fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        fout.write("\t".join(vcf_info.out_samples) + '\n')

        for entry_set in monitored_pool(to_call, jobs, threads):
            fout.write(entry_set)

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
    parser.add_argument("input", metavar="INPUT", type=str, nargs="+",
                        help="Input VCFs to harmonize")
    parser.add_argument("-r", "--region", type=str, required=True,
                        help="Bed filename or comma-separated list of chrom:start-end regions to process")
    parser.add_argument("-f", "--reference", type=str, required=True,
                        help="Reference")
    parser.add_argument("-o", "--output", type=str, default="phab_out.vcf.gz",
                        help="Output VCF")
    parser.add_argument("--maxsize", type=int, default=50000,
                        help="Maximum size of variant to incorporate into haplotypes (%(default)s)")
    parser.add_argument("--passonly", action="store_true",
                        help="Only incorporate passing variants (%(default)s)")
    parser.add_argument("--dedup", action="store_true",
                        help="Dedupicate haplotypes before MSA")
    parser.add_argument("--buffer", type=int, default=100,
                        help="Number of reference bases before/after region to add to MSA (%(default)s)")
    parser.add_argument("--align", type=str, choices=["mafft", "wfa", "poa"], default="poa",
                        help="Alignment method (%(default)s)")
    parser.add_argument("--mafft-params", type=str, default=DEFAULT_MAFFT_PARAM,
                        help="Parameters for mafft, wrap in a single quote (%(default)s)")
    parser.add_argument("--samples", type=str, default=None,
                        help="Subset of samples to MSA")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of threads (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    if args.samples is not None:
        args.samples = args.samples.split(',')

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

    if not os.path.exists(args.reference):
        logging.error("Reference %s does not exist", args.reference)
        check_fail = True

    for i in args.input:
        if not os.path.exists(i):
            logging.error("File %s does not exist", i)
            check_fail = True
        else:
            if not i.endswith(".gz"):
                logging.error("File %s does not end with .gz. Must be bgzip'd", i)
                check_fail = True
            if not truvari.check_vcf_index(i):
                logging.error("File %s must be indexed", i)
                check_fail = True
    if not args.samples is None and not len(args.samples) == len(set(args.samples)):
        logging.error("Redundant sample names in --samples")
        check_fail = True
    return check_fail


def get_align_method(method, params=DEFAULT_MAFFT_PARAM):
    """
    Helper for organizing how to get the different alignment methods
    """
    if method == "mafft":
        align_method = functools.partial(run_mafft, params=params)
    elif method == "wfa":
        align_method = run_wfa
    elif method == "poa":
        align_method = run_poa
    else:
        raise ValueError(f"Unknown alignment method {method}")
    return align_method


def phab_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)
    truvari.setup_logging(args.debug, show_version=True)
    if check_requirements(args.align) or check_params(args):
        logging.error("Couldn't run Truvari. Please fix parameters\n")
        sys.exit(100)

    m_vcf_info = VCFtoHaplotypes(args.reference,
                                 args.input, args.samples,
                                 passonly=args.passonly,
                                 max_size=args.maxsize)

    all_regions = parse_regions(args.region)
    method = get_align_method(args.align, args.mafft_params)

    run_phab(m_vcf_info, all_regions, args.output,
             buffer=args.buffer,
             align_method=method,
             dedup=args.dedup,
             threads=args.threads)

    logging.info("Finished phab")
