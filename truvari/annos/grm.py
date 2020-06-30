"""
Maps graph edge kmers with BWA to assess Graph Reference Mappability
"""
import os
import re
import queue
import logging
import argparse
import itertools

import pysam
import joblib
import truvari
import pandas as pd

try:
    from bwapy import BwaAligner
    HASBWALIB=True
except ModuleNotFoundError:
    HASBWALIB=False

def index(vcf, bufcnt=1, chunksize=100000):
    """
    Get the positions of chunks, yield them as they're created
    """
    v = pysam.VariantFile(vcf, threads=2)
    que = queue.Queue()
    cur = []
    cnt = 0
    for entry in v:
        cur.append(entry)
        if len(cur) >= chunksize:
            que.put(cur)
            cnt += 1
            cur = []
            if que.qsize() >= bufcnt:
                yield que.get()
    cnt += 1
    que.put(cur)
    while not que.empty():
        yield que.get()
    logging.info("%d chunks", cnt)


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
        up = ref.fetch(entry.chrom, start - kmer, start + kmer)
        dn = ref.fetch(entry.chrom, end - kmer, end + kmer)
    except Exception:
        return None
    seq = entry.alts[0]
    if seq[0] == up[kmer]:
        seq = seq[1:]  # trimming...
    # alternate
    hap = up[:kmer] + seq + dn[-kmer:]

    return up, dn, hap[:kmer * 2], hap[-kmer * 2:]


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
    Maps the kmer and returns the 
    max/min
    if chrom/pos is provided:
        remove any hits that maps over this position. This is a filter of the reference - we want to know where *else* it hits
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
    parser = argparse.ArgumentParser(prog="kmap", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", required=True,
                        help="Input VCF")
    parser.add_argument("-r", "--reference", required=True,
                        help="BWA indexed reference")
    parser.add_argument("-o", "--output", default="results.jl",
                        help="Output dataframe (%(default)s)")
    parser.add_argument("-k", "--kmersize", default=50, type=int,
                        help="Size of kmer to map (%(default)s)")
    parser.add_argument("-t", "--threads", default=os.cpu_count(), type=int,
                        help="Number of threads (%(default)s)")
    args = parser.parse_args(args)
    return args


def process_entry(entries, ref, aligner, kmersize):
    rows = []
    for entry in entries:
        if "SVLEN" not in entry.info:
            continue
        kmers = make_kmers(ref, entry, kmersize//2)
        if kmers is None:
            continue
        ref_up, ref_dn, alt_up, alt_dn = kmers
        ty = truvari.entry_variant_type(entry)
        sz = truvari.entry_size(entry)
        rcov, acov = entry.samples[0]["AD"]
        try:
            gref, ghet, ghom = entry.samples[0]["PL"]
        except Exception:
            gref, ghet, ghom = 0, 0, 0

        result = ["%s:%d-%d.%s" % (entry.chrom, entry.start, entry.stop, entry.alts[0]),
                  ty, sz, rcov, acov, gref, ghet, ghom]

        if ty == "INS":
            # Only want a single reference kmer
            ref_stats = map_stats(aligner, ref_up, entry.chrom, entry.start)
            result.extend(ref_stats + ref_stats)
            result.extend(map_stats(aligner, alt_up))
            result.extend(map_stats(aligner, alt_dn))
        elif ty == "DEL":
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

    return rows


def grm_main(cmdargs):
    """
    Builds a graph-genome from the vcf and reference,
    creates the sets of kmers at edges in the graph
    maps those kmers globally to the reference
    reports mapping metrics

    Todo: 
    - document the bwa package and how to install
    - better column names along with documentation

    Long-term: 
    This should probably work to annotate the VCF's INFO fields
    """
    if not HASBWALIB:
        logging.error("bwapy not available. Please install https://github.com/nanoporetech/bwapy")
        exit(1)

    args = parse_args(cmdargs)
    ref = pysam.FastaFile(args.reference)
    aligner = BwaAligner(args.reference)

    p_run = joblib.Parallel(n_jobs=args.threads, prefer="threads")

    logging.info("Indexing")
    events = index(args.input, args.threads)
    logging.info("Processing")
    chunks = p_run(joblib.delayed(process_entry)(e, ref, aligner, args.kmersize) for e in events)
    header = ["key", "svtype", "svlen", "ref_cov", "alt_cov", "gref", "ghet", "ghom"]
    for prefix in ["rup_", "rdn_", "aup_", "adn_"]:
        for key in ["nhits", "avg_q", "avg_ed", "avg_mat", "avg_mis", "dir_hits", "com_hits", "max_q",
                    "max_ed", "max_mat", "max_mis", "max_strand",
                    "min_q", "min_ed", "min_mat", "min_mis", "min_strand"]:
            header.append(prefix + key)

    logging.info("Saving")
    data = pd.DataFrame(itertools.chain(*chunks), columns=header)
    joblib.dump(data, args.output)
    logging.info("df shape %s", data.shape)
    logging.info("Finished")
