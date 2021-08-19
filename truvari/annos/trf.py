"""
Intersect vcf with reference simple repeats and report
how the an alternate allele affects the copies using TRF
"""
import os
import sys
import types
import decimal
import logging
import argparse
import tempfile
from collections import defaultdict

import pysam
from acebinf import cmd_exe, setup_logging
import truvari

trfshared = types.SimpleNamespace()

def parse_args(args):
    """
    Pull the command line parameters
    """
    #def restricted_float(x):
    #    x = float(x)
    #    if x < 0.0 or x > 1.0:
    #        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    #    return x
    parser = argparse.ArgumentParser(prog="trf", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, default="/dev/stdin",
                        help="VCF to annotate (stdin)")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output filename (stdout)")
    parser.add_argument("-e", "--executable", type=str, default="trf409.linux64",
                        help="Path to tandem repeat finder (%(default)s)")
    parser.add_argument("-T", "--trf-params", type=str, default="2 7 7 80 10 50 500 -m -f -h -d -ngs",
                        help="Default parameters to send to trf (%(default)s)")
    parser.add_argument("-s", "--simple-repeats", type=str,
                        help="Simple repeats bed")
    parser.add_argument("-f", "--reference", type=str,
                        help="Reference fasta file")
    #parser.add_argument("-m", "--min-length", type=int, default=50,
                        #help="Minimum size of entry to annotate (%(default)s)")
    #parser.add_argument("-t", "--threshold", type=restricted_float, default=.8,
                        #help="Threshold for pct of allele covered (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    setup_logging(args.debug)
    return args

class TRFAnno():
    """ Class for trf annotation """
    def __init__(self, executable="trf409.linux64",
                 trf_params="2 7 7 80 10 50 500 -m -f -h -d -ngs",
                 tmpdir=None, simple_repeats=None, reference=None):
        """ setup """
        # This is slow, so I only want to do it once
        # And actually, let's put it in a shared memory spot just
        # because we can
        self.simple_repeats = simple_repeats if simple_repeats else trfshared.simple_repeats
        self.reference = pysam.FastaFile(reference) if reference else pysam.FastaFile(trfshared.reference)
        self.executable = executable
        if "-ngs" not in trf_params:
            trf_params = trf_params + " -ngs "
        self.trf_params = trf_params
        if tmpdir is None:
            tmpdir = tempfile._get_default_tempdir()
        # Where we write the fasta entries
        self.fa_fn = os.path.join(tmpdir, next(tempfile._get_candidate_names()))
        self.tr_fn = os.path.join(tmpdir, next(tempfile._get_candidate_names()))
        # lookup from the vcf entries to the hits
        self.srep_lookup = defaultdict(list)

    def make_seq(self, srep, entry):
        """
        Make the haplotype sequence
        """
        start = srep.begin
        end = srep.end
        if entry.info["SVTYPE"] == "INS":
            ref_seq = self.reference.fetch(entry.chrom, start, end)
            m_seq = ref_seq[:entry.start - start] + entry.alts[0] + ref_seq[entry.stop - start:]
        elif entry.info["SVTYPE"] == "DEL":
            # Only want to consider the srep region
            m_start = max(start, entry.start)
            m_end = min(end, entry.stop)
            m_seq = self.reference.fetch(entry.chrom, m_start, m_end)
        else:
            logging.critical("Can only consider entries with 'SVTYPE' INS/DEL")
            sys.exit(1)
        return m_seq

    def calc_copy_diff(self, key, trf_hits):
        """
        trf_hits are on the alternate
        srep_lookup is on the reference

        Given the variant's key and the data from the trf hits
        pick the best one
        data is a dictionary of the hits with a key of the hits' sequence
        seqs is the reference simple repeat hits
        """
        found_len = 0
        copy_diff = None
        logging.debug('variant: %s', key)
        logging.debug('trf_hits:')
        for i in trf_hits.values():
            logging.debug(i)
        logging.debug('srep_hits:')
        for srep in self.srep_lookup[key]:
            logging.debug(srep)
            repeat = srep["SREP_repeats"][0]
            copies = srep["SREP_copies"][0]
            if repeat in sorted(trf_hits.keys(), reverse=True):
                if len(repeat) >= found_len:
                    copy_diff = trf_hits[repeat]['copies'] - copies
        if copy_diff:
            logging.debug('cdiff %.1f', copy_diff)
        else:
            logging.debug('no hit')
        return copy_diff

    def run_trf(self, entries):
        """
        Given a list of entries, run TRF on each of them
        """
        # Write seqs
        with open(self.fa_fn, 'w') as fout:
            for entry in entries:
                hits = self.simple_repeats[entry.chrom].overlap(entry.start, entry.stop)
                if not hits:
                    continue
                key = f"{entry.chrom}:{entry.start}-{entry.stop}.{hash(entry.alts[0])}"
                for srep in hits:
                    self.srep_lookup[key].append(srep.data)
                    seq = self.make_seq(srep, entry)
                fout.write(f">{key}\n{seq}\n")

        # Run it
        ret = cmd_exe(f"{self.executable} {self.fa_fn} {self.trf_params} > {self.tr_fn}")
        if ret.ret_code != 0:
            logging.error("Couldn't run trf")
            logging.error(str(ret))
            sys.exit(ret.ret_code)

        #Then I need to parse and return
        return self.parse_output()

    def parse_output(self):
        """
        Parse the outputs from trf, turn to a dictionary
        """
        trf_cols = [("start", int),
                    ("end", int),
                    ("period", int),
                    ("copies", float),
                    ("consize", int),
                    ("pctmat", int),
                    ("pctindel", int),
                    ("score", int),
                    ("A", int),
                    ("C", int),
                    ("G", int),
                    ("T",  int),
                    ("entropy", float),
                    ("repeat", str),
                    ("unk1", None),
                    ("unk2", None),
                    ("unk3", None)]
        hits = {}
        with open(self.tr_fn, 'r') as fh:
            name = fh.readline()
            if name == "": # no hits
                return hits
            name = name.strip()[1:]
            trf_hits = {} # all the seqence hits we get for this entry
            while True:
                line = fh.readline()
                if line == "":
                    break
                if line.startswith("@"):
                    key = name.strip()
                    copy_diff = self.calc_copy_diff(key, trf_hits)
                    if copy_diff:
                        hits[key] = copy_diff
                    name = line.strip()[1:]
                    continue
                line = line.strip().split(' ')
                data = {x[0]: x[1](y) for x, y in zip(trf_cols, line) if not x[0].startswith("unk")}
                trf_hits[data["repeat"]] = data
        return hits

def edit_header(header):
    """
    New VCF INFO fields
    """
    header = header.copy()
    # if intersect_only: Do I want to sometimes turn this off?
    # Probably, actully
    header.add_line(('##INFO=<ID=SimpleRepeatDiff,Number=1,Type=Float,'
                     'Description="TRF simple repeat copy difference">'))
    return header

def trf_main(cmdargs):
    """ TRF annotation """
    args = parse_args(cmdargs)
    trfshared.reference = args.reference
    #"/home/english/truvari/repo_utils/test_files/reference.fa"
    #refanno = "/home/english/truvari/repo_utils/test_files/simplerepeat.txt.gz"
    tree = truvari.make_bedanno_tree(args.simple_repeats)
    trfshared.simple_repeats = tree[0]
    #vcf = "/home/english/truvari/repo_utils/test_files/multi.vcf.gz"
    v = pysam.VariantFile(args.input)
    to_consider = []
    for entry in v:
        # need to set a minmum length probably
        if "SVTYPE" in entry.info: #and entry.info["SVTYPE"] == "INS":
            to_consider.append(entry)
    tanno = TRFAnno(executable=args.executable)
    #/home/english/truvari/repo_utils/test_files/external/trf
    annos = tanno.run_trf(to_consider)

    v = pysam.VariantFile(args.input)
    new_header = edit_header(v.header)
    o = pysam.VariantFile(args.output, 'w', header=new_header)
    decimal.getcontext().prec = 1
    for entry in v:
        if "SVTYPE" in entry.info:
            key = f"{entry.chrom}:{entry.start}-{entry.stop}.{hash(entry.alts[0])}"
            if key in annos:
                entry = truvari.copy_entry(entry, new_header)
                entry.info["SimpleRepeatDiff"] = annos[key]
        o.write(entry)
    logging.info("Finished trf")

if __name__ == '__main__':
    trf_main(sys.argv[1:])
