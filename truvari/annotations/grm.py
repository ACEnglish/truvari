"""
Maps graph edge kmers with BWA to assess Graph Reference Mappability
"""
# pylint: disable=too-many-locals
import sys
import re
import types
import logging
import argparse
import functools
import multiprocessing
from collections import namedtuple

import pysam
import tabix
import joblib
import pandas as pd

try:
    from bwapy import BwaAligner
    HASBWALIB = True
except (OSError,  ModuleNotFoundError):
    HASBWALIB = False

import truvari

try:
    from setproctitle import setproctitle  # pylint: disable=import-error,useless-suppression
except ModuleNotFoundError:
    def setproctitle(_):
        """ dummy function """
        return


def make_kmers(ref, entry, kmer=25):
    """
    Make ref/alt kmers
    kmer should be half the kmer size you're actually trying to make
    Returns 4 kmers
    ref up, ref dn, alt up, alt dn
    """
    start = entry.start
    end = entry.stop
    try:
        reflen = ref.get_reference_length(entry.chrom)
        up = ref.fetch(entry.chrom, max(start - kmer, 0),
                       min(start + kmer, reflen))
        dn = ref.fetch(entry.chrom, max(end - kmer, 0),
                       min(end + kmer, reflen))

        if len(up) < kmer or len(dn) < kmer:
            return None
        seq = entry.alts[0]
        if len(seq) >= 1 and seq[0] == up[kmer]:
            seq = seq[1:]  # trimming...
        # alternate
        hap = up[:kmer] + seq + dn[-kmer:]
        return up, dn, hap[:kmer * 2], hap[-kmer * 2:]
    except Exception as e:  # pylint: disable=broad-except
        logging.warning(f"{e} for {str(entry)[:20]}...")
        return None


cigmatch = re.compile("[0-9]+[MIDNSHPX=]")


def cig_pctsim(cigar):
    """
    Cheap way to count how many non-matches there are.
    Note M just means the bases are Matched up in the alignment,
    not that the matched bases are the same. Therfore, M's could
    be "mis-matches".
    """
    match = 0
    miss = 0
    for i in cigmatch.findall(cigar):
        if i[-1] == "M":
            match += int(i[:-1])
        else:
            miss += int(i[:-1])
    return match, miss


def map_stats(aligner, kmer, chrom=None, pos=None):
    """
    Maps the kmer and returns the max/min
    if chrom/pos is provided remove any hits that maps over this position. This
    is a filter of the reference - we want to know where *else* it hits
    """
    nhits = 0
    dir_hits = 0
    com_hits = 0
    max_q = -1
    max_strand = 0
    max_ed = 0
    max_mat = 0
    max_mis = 0

    min_q = 1000
    min_strand = 0
    min_ed = 0
    min_mat = 0
    min_mis = 0

    avg_q = 0
    avg_ed = 0
    avg_mat = 0
    avg_mis = 0

    for aln in aligner.align_seq(kmer):
        if chrom is not None and aln.rname == chrom and aln.pos < pos < aln.pos + len(kmer):
            continue  # Skipping reference hitting where it came from
        match, miss = cig_pctsim(aln.cigar)

        nhits += 1
        avg_q += aln.mapq
        avg_ed += aln.NM
        avg_mat += match
        avg_mis += miss

        if aln.orient == "+":
            dir_hits += 1
        else:
            com_hits += 1

        # MAX
        if aln.mapq > max_q:
            max_q = aln.mapq
            max_strand = aln.orient
            max_ed = aln.NM
            max_mat = match
            max_mis = miss

        # MIN
        if aln.mapq < min_q:
            min_q = aln.mapq
            min_strand = aln.orient
            min_ed = aln.NM
            min_mat = match
            min_mis = miss

    if nhits != 0:
        avg_q /= nhits
        avg_ed /= nhits
        avg_mat /= nhits
        avg_mis /= nhits
    return nhits, avg_q, avg_ed, avg_mat, avg_mis, dir_hits, com_hits, \
        max_q, max_ed, max_mat, max_mis, max_strand, \
        min_q, min_ed, min_mat, min_mis, min_strand


def parse_args(args):
    """
    Argument parsing
    """
    parser = argparse.ArgumentParser(prog="grm", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", required=True,
                        help="Input VCF")
    parser.add_argument("-r", "--reference", required=True,
                        help="BWA indexed reference")
    parser.add_argument("-R", "--regions", default=None,
                        help="Bed file of regions to parse (None)")
    parser.add_argument("-o", "--output", default="results.jl",
                        help="Output dataframe (%(default)s)")
    parser.add_argument("-k", "--kmersize", default=50, type=truvari.restricted_int,
                        help="Size of kmer to map (%(default)s)")
    parser.add_argument("-m", "--min-size", default=25, type=truvari.restricted_int,
                        help="Minimum size of variants to map (%(default)s)")
    parser.add_argument("-t", "--threads", default=1, type=truvari.restricted_int,
                        help="Number of threads (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug, show_version=True)
    return args


def parse_infos(infos):
    """
    Convert INFO field to dict
    """
    for info in infos:
        s = info.split("=", 1)
        if len(s) == 1:
            yield info, None
        else:
            yield s


Entry = namedtuple("Entry", ["chrom", "start", "stop", "ref", "alts", "info"])


def line_to_entry(fields):
    """
    Replicate pysam entries, but lighter
    """
    if len(fields) < 8:
        raise RuntimeError(f"Can't parse line: {fields}")
    info_field = fields[7].split(";")
    info_dict = dict(parse_infos(info_field))
    ref = fields[3]
    alts = fields[4].split(",")
    start = int(fields[1]) - 1
    return truvari.VariantRecord(Entry(fields[0],  # chrom
                 start,
                 start + len(ref),  # stop
                 ref,
                 alts,
                 info_dict))


def read_vcf_lines(in_fn, ref_name, start, stop):
    """
    Faster VCF parsing
    """
    logging.debug(f"Starting region {ref_name}:{start}-{stop}")

    tb = tabix.open(in_fn)
    try:
        yield from tb.query(ref_name, start, stop)
    except tabix.TabixError as e:
        logging.warning(f"Region {ref_name}:{start}-{stop} failed: {e}")
    setproctitle(f"grm done {ref_name}:{start}-{stop}")
    logging.debug(f"Done region {ref_name}:{start}-{stop}")


def process_entries(ref_section, grm_shared):
    """
    Calculate GRMs for a set of vcf entries
    """
    grm_shared.aligner = BwaAligner(grm_shared.ref_filename, '-a')
    ref_name, start, stop = ref_section
    ref = pysam.FastaFile(grm_shared.ref_filename)
    aligner = grm_shared.aligner
    kmersize = grm_shared.kmersize
    header = grm_shared.header
    minsize = grm_shared.min_size
    rows = []
    next_progress = 0
    for line in read_vcf_lines(grm_shared.input, ref_name, start, stop):
        if "SVLEN" not in line[7] and abs(len(line[3]) - len(line[4])) < minsize:
            continue
        entry = line_to_entry(line)
        if next_progress == 0:
            next_progress = 1000
            setproctitle(
                f"grm processing {entry.chrom}:{entry.start} in {ref_name}:{start}-{stop}")
        else:
            next_progress -= 1
        if "SVLEN" not in entry.info:
            continue

        kmers = make_kmers(ref, entry, kmersize//2)
        if kmers is None:
            continue

        ref_up, ref_dn, alt_up, alt_dn = kmers
        ty = entry.var_type()

        result = ["%s:%d-%d.%s" %  # pylint: disable=consider-using-f-string
                  (entry.chrom, entry.start, entry.stop, entry.alts[0])]
        if ty == truvari.SV.INS:
            # Only want a single reference kmer
            ref_stats = map_stats(aligner, ref_up, entry.chrom, entry.start)
            result.extend(ref_stats + ref_stats)
            result.extend(map_stats(aligner, alt_up))
            result.extend(map_stats(aligner, alt_dn))
        elif ty == truvari.SV.DEL:
            result.extend(map_stats(aligner, ref_up, entry.chrom, entry.start))
            result.extend(map_stats(aligner, ref_dn, entry.chrom, entry.stop))
            # Only want a single alternate kmer
            alt_stats = map_stats(aligner, alt_up)
            result.extend(alt_stats + alt_stats)
        else:
            result.extend(map_stats(aligner, ref_up, entry.chrom, entry.start))
            result.extend(map_stats(aligner, ref_dn, entry.chrom, entry.stop))
            result.extend(map_stats(aligner, alt_up))
            result.extend(map_stats(aligner, alt_dn))
        rows.append(result)
    data = pd.DataFrame(rows, columns=header)
    logging.debug(f"Chunk {ref_section} finished, shape={data.shape}")
    return data


def grm_main(cmdargs):
    """
    Builds a graph-genome from the vcf and reference,
    creates the sets of kmers at edges in the graph
    maps those kmers globally to the reference
    reports mapping metrics

    Todo:
    - document the bwa package and how to install
    - better column names along with documentation
    """
    if not HASBWALIB:
        logging.error("bwapy isn't available on this machine")
        sys.exit(1)
    args = parse_args(cmdargs)
    if not args.regions:
        m_ranges = truvari.ref_ranges(args.reference)
    else:
        m_ranges = truvari.bed_ranges(args.regions)

    header = ["key"]
    for prefix in ["rup_", "rdn_", "aup_", "adn_"]:
        for key in ["nhits", "avg_q", "avg_ed", "avg_mat", "avg_mis", "dir_hits", "com_hits", "max_q",
                    "max_ed", "max_mat", "max_mis", "max_strand",
                    "min_q", "min_ed", "min_mat", "min_mis", "min_strand"]:
            header.append(prefix + key)
    grm_shared = types.SimpleNamespace()
    grm_shared.header = header
    grm_shared.ref_filename = args.reference
    grm_shared.kmersize = args.kmersize
    grm_shared.input = args.input
    grm_shared.min_size = args.min_size
    m_process_entries = functools.partial(
        process_entries, grm_shared=grm_shared)
    with multiprocessing.Pool(args.threads, maxtasksperchild=1) as pool:
        logging.info("Processing")
        chunks = pool.imap(m_process_entries, m_ranges)
        pool.close()
        data = pd.concat(chunks, ignore_index=True)
        logging.info("Saving; df shape %s", data.shape)
        joblib.dump(data, args.output)
        logging.info("Finished grm")
        pool.join()
