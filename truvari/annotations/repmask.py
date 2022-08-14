""" Wrapper around RepeatMasker to annotate insertion sequences in a VCF """
import os
import sys
import logging
import argparse
import tempfile
from collections import defaultdict

import pysam
import truvari

DEFAULTPARAMS = "-pa {threads} -qq -e hmmer -species human -lcambig -nocut -div 50 -no_id -s {fasta}"


def paren_int(number):
    """
    returns an integer from a string "([0-9]+)"
    """
    return int(number.replace('(', '').replace(')', ''))


class RepMask():
    """ Class for RepeatMasker annotation """
    REPCOLS = [("RM_score", int),
               ("RM_pdiv", float),
               ("RM_pdel", float),
               ("RM_pins", float),
               ("RM_query", str),
               ("RM_qstart", int),
               ("RM_qend", int),
               ("RM_qleft", paren_int),  # need custom here (it's "(int)")
               ("RM_strand", str),
               ("RM_repeat", str),
               ("RM_clsfam", str),  # need custom here.. usually class/family
               ("RM_tstart", paren_int),
               ("RM_tend", paren_int),
               ("RM_tleft", str)]

    def __init__(self, in_vcf, out_vcf="/dev/stdout", executable="RepeatMasker",
                 min_length=50, max_length=50000, threshold=0.8,  rm_params=DEFAULTPARAMS, threads=1):
        """ The setup """
        self.in_vcf = str(in_vcf)
        self.out_vcf = out_vcf
        self.executable = executable
        self.min_length = min_length
        self.max_length = max_length
        self.threshold = threshold
        self.rm_params = rm_params
        self.threads = threads
        self.n_header = None
        self.cmd = f"{self.executable} {self.rm_params}"

    def edit_header(self, header):
        """
        Edits and holds on to the header
        """
        header.add_line(('##INFO=<ID=RM_score,Number=1,Type=Integer,'
                         'Description="RepMask bit score">'))
        header.add_line(('##INFO=<ID=RM_repeat,Number=1,Type=String,'
                         'Description="RepMask matching repeat">'))
        header.add_line(('##INFO=<ID=RM_clsfam,Number=1,Type=String,'
                         'Description="RepMask repeat class/family ">'))
        # TODO: Need to put a source line that says this thing was run with whatever parameters
        self.n_header = header

    def extract_seqs(self):
        """
        Create the fasta file of all the sequences
        Returns the fasta file name
        """
        ret = tempfile.NamedTemporaryFile(mode='w', delete=False) # pylint: disable=consider-using-with
        tot_cnt = 0
        cnt = 0
        cntbp = 0
        with pysam.VariantFile(self.in_vcf) as fh:
            for pos, entry in enumerate(fh):
                tot_cnt += 1
                entry_size = truvari.entry_size(entry)
                if self.min_length <= entry_size <= self.max_length:
                    cnt += 1
                    cntbp += entry_size
                    if truvari.entry_variant_type(entry) == "INS":
                        ret.write(f">{pos}\n{entry.alts[0]}\n")
                    else:
                        ret.write(f">{pos}\n{entry.ref}\n")
        logging.info(
            f"Extracted {cnt} sequences ({cntbp}bp) from {tot_cnt} entries")
        return ret.name

    @staticmethod
    def parse_output(faout):
        """
        Parses the RepeatMasker output
        """
        hits = defaultdict(list)
        with open(faout, 'r') as fh:
            # header lines
            fh.readline()
            fh.readline()
            fh.readline()
            for line in fh:
                data = line.strip().split()
                data = {x[0]: x[1](y) for x, y in zip(RepMask.REPCOLS, data)}
                hits[data["RM_query"]].append(data)
        return hits

    def annotate_seqs(self, fasta):
        """
        Runs repeat masker on a fasta
        Parses the {fasta}.out
        Returns dict of {entry_pos_id: [{hitdict}, ...], ...}
        """
        logging.info("Starting RepeatMasker")
        cmd = self.cmd.format(threads=self.threads, fasta=fasta)
        ret = truvari.cmd_exe(cmd)
        if ret.ret_code != 0:
            logging.error("Couldn't run RepeatMasker")
            logging.error(str(ret))
            sys.exit(ret.ret_code)
        logging.info("Finished RepeatMasker")

        hits = self.parse_output(f"{fasta}.out")
        return hits

    def annotate_entry(self, entry, hits):
        """
        Annotates a single entry with given hits
        """
        best_hit_pct = 0
        best_hit = None
        entry_size = truvari.entry_size(entry)
        for hit in hits:
            size_aln = abs(hit["RM_qstart"] - hit["RM_qend"]) + 1
            pct = size_aln / entry_size  # The TR that covers the most of the sequence
            # I'm taking the single best... So I might be 'under annotating'
            # Also, I might want to consider the score?
            if pct >= self.threshold and pct > best_hit_pct:
                best_hit_pct = pct
                best_hit = hit
        return self.edit_entry(entry, best_hit)

    def annotate_vcf(self):
        """
        Annotates entries in the vcf and writes to new vcf
        """
        hits = self.annotate_seqs(self.extract_seqs())
        with pysam.VariantFile(self.in_vcf, 'r') as fh:
            self.edit_header(fh.header.copy())
            with pysam.VariantFile(self.out_vcf, 'w', header=self.n_header) as out:
                for pos, entry in enumerate(fh):
                    pos = str(pos)
                    if pos in hits:
                        entry = self.annotate_entry(entry, hits[pos])
                    out.write(entry)

    def edit_entry(self, entry, rm_hit):
        """
        puts the annos in vcf entry
        return the edited entry
        """
        if not rm_hit:
            return entry
        entry.translate(self.n_header)

        entry.info["RM_score"] = rm_hit["RM_score"]
        entry.info["RM_repeat"] = rm_hit["RM_repeat"]
        entry.info["RM_clsfam"] = rm_hit["RM_clsfam"]

        return entry


def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="repmask", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", type=str, required=True,
                        help="VCF to annotate (%(default)s)")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output filename (%(default)s)")
    parser.add_argument("-e", "--executable", type=str, default="RepeatMasker",
                        help="Path to RepeatMasker (%(default)s)")
    parser.add_argument("-m", "--min-length", type=truvari.restricted_int, default=50,
                        help="Minimum size of entry to annotate (%(default)s)")
    parser.add_argument("-M", "--max-length", type=truvari.restricted_int, default=50000,
                        help="Maximum size of entry to annotate (%(default)s)")
    parser.add_argument("-t", "--threshold", type=truvari.restricted_float, default=.8,
                        help="Threshold for pct of allele covered (%(default)s)")
    parser.add_argument("-p", "--params", type=str, default=DEFAULTPARAMS,
                        help="Default parameter string to send to RepeatMasker (%(default)s)")
    parser.add_argument("-T", "--threads", type=truvari.restricted_int, default=os.cpu_count(),
                        help="Number of threads to use (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    return args


def rmk_main(cmdargs):
    """ Main """
    args = parse_args(cmdargs)
    anno = RepMask(in_vcf=args.input,
                   out_vcf=args.output,
                   executable=args.executable,
                   min_length=args.min_length,
                   max_length=args.max_length,
                   threshold=args.threshold,
                   rm_params=args.params,
                   threads=args.threads)
    anno.annotate_vcf()
    logging.info("Finished repmask")
