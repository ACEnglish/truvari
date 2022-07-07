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
import multiprocessing
from io import StringIO
from collections import defaultdict

import pysam
import tabix
import truvari

trfshared = types.SimpleNamespace()

try:
    from setproctitle import setproctitle  # pylint: disable=import-error,useless-suppression
except ModuleNotFoundError:
    def setproctitle(_):
        """ dummy function """
        return


def fetch_simple_repeats(chrom, start, stop):
    """
    Parse a simple repeats bed file and yield them
    """
    header = [('chrom', str), ('start', int), ('end', int), ('period', float),
              ('copies', float), ('score', int), ('entropy', float), ('repeat', str)]
    tb = tabix.open(trfshared.args.simple_repeats)
    try:
        for i in tb.query(chrom, start, stop):
            yield {fmt[0]: fmt[1](x) for x, fmt in zip(i, header)}
    except tabix.TabixError as e:
        logging.warning(f"Region {chrom}:{start}-{stop} failed: {e}")
    return False


def process_entries(ref_section):
    """
    Process vcf lines from a reference section
    """
    chrom, start, stop = ref_section
    logging.debug(f"Starting region {chrom}:{start}-{stop}")
    setproctitle(f"trf {chrom}:{start}-{stop}")
    vcf = pysam.VariantFile(trfshared.args.input)

    to_consider = []
    for entry in vcf.fetch(chrom, start, stop):
        # Prevent duplication
        if not (entry.start >= start and entry.start < stop):
            continue
        if truvari.entry_size(entry) >= trfshared.args.min_length:
            to_consider.append(entry)

    if not to_consider:
        return (chrom, start, stop, "")

    tanno = TRFAnno(executable=trfshared.args.executable,
                    trf_params=trfshared.args.trf_params)
    tanno.run_trf(to_consider)

    v = pysam.VariantFile(trfshared.args.input)
    new_header = edit_header(v.header)
    out = StringIO()
    decimal.getcontext().prec = 1
    for entry in v.fetch(chrom, start, stop):
        # Prevent duplication
        if not (entry.start >= start and entry.start < stop):
            continue
        if truvari.entry_size(entry) >= trfshared.args.min_length:
            key = f"{entry.chrom}:{entry.start}-{entry.stop}.{hash(entry.alts[0])}"
            entry = tanno.annotate(entry, key, new_header)
        out.write(str(entry))
    out.seek(0)
    setproctitle(f"trf done {chrom}:{start}-{stop}")
    logging.debug(f"Done region {chrom}:{start}-{stop}")
    return (chrom, start, stop, out.read())


def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="trf", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="VCF to annotate")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output filename (stdout)")
    parser.add_argument("-e", "--executable", type=str, default="trf409.linux64",
                        help="Path to tandem repeat finder (%(default)s)")
    parser.add_argument("-T", "--trf-params", type=str, default="3 7 7 80 5 40 500 -h -ngs",
                        help="Default parameters to send to trf (%(default)s)")
    parser.add_argument("-s", "--simple-repeats", type=str, required=True,
                        help="Simple repeats bed")
    parser.add_argument("-f", "--reference", type=str, required=True,
                        help="Reference fasta file")
    parser.add_argument("-m", "--min-length", type=truvari.restricted_int, default=50,
                        help="Minimum size of entry to annotate (%(default)s)")
    parser.add_argument("-M", "--max-length", type=truvari.restricted_int, default=10000,
                        help="Maximum size of sequence to run through trf (%(default)s)")
    parser.add_argument("-t", "--threads", type=truvari.restricted_int, default=multiprocessing.cpu_count(),
                        help="Number of threads to use (%(default)s)")
    parser.add_argument("-C", "--chunk-size", type=truvari.restricted_int, default=1,
                        help="Size (in mbs) of reference chunks for parallelization (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")

    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
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
        self.simple_repeats = simple_repeats if simple_repeats else trfshared.args.simple_repeats
        self.reference = pysam.FastaFile(
            reference) if reference else pysam.FastaFile(trfshared.args.reference)
        self.executable = executable
        if "-ngs" not in trf_params:
            trf_params = trf_params + " -ngs "
        self.trf_params = trf_params
        if tmpdir is None:
            tmpdir = tempfile._get_default_tempdir()
        # Where we write the fasta entries
        self.fa_fn = os.path.join(tmpdir, next(
            tempfile._get_candidate_names()))
        self.tr_fn = os.path.join(tmpdir, next(
            tempfile._get_candidate_names()))
        # lookup from the vcf entries to the hits
        self.trf_lookup = defaultdict(dict)
        self.srep_lookup = defaultdict(dict)

    def make_seq(self, start, end, entry):
        """
        Make the haplotype sequence
        """
        svtype = truvari.entry_variant_type(entry)
        if svtype == "INS":
            ref_seq = self.reference.fetch(entry.chrom, start, end)
            m_seq = ref_seq[:entry.start - start] + \
                entry.alts[0] + ref_seq[entry.stop - start:]
        elif svtype == "DEL":
            m_start = max(start, entry.start)
            m_end = min(end, entry.stop)
            m_seq = self.reference.fetch(entry.chrom, m_start, m_end)
        else:
            logging.critical("Can only consider entries with 'SVTYPE' INS/DEL")
            sys.exit(1)
        return m_seq

    def run_trf(self, entries):
        """
        Given a list of entries, run TRF on each of them
        """
        n_seqs = 0
        with open(self.fa_fn, 'w') as fout:
            for entry in entries:
                hits = list(fetch_simple_repeats(
                    entry.chrom, entry.start - 1, entry.stop + 1))
                if not hits:
                    continue
                key = f"{entry.chrom}:{entry.start}-{entry.stop}.{hash(entry.alts[0])}"
                for srep in hits:
                    seq = self.make_seq(srep["start"], srep["end"], entry)
                    self.srep_lookup[key][srep['repeat']] = srep
                    if len(seq) <= trfshared.args.max_length:
                        n_seqs += 1
                        fout.write(f">{key}\n{seq}\n")

        if not n_seqs:
            return

        # Run it
        ret = truvari.cmd_exe(
            f"{self.executable} {self.fa_fn} {self.trf_params} > {self.tr_fn}")
        if ret.ret_code != 0:
            logging.error("Couldn't run trf. Check Parameters")
            logging.error(f"{self.executable} {self.fa_fn} {self.trf_params} > {self.tr_fn}")
            logging.error(str(ret))
            return

        self.parse_output()

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
        with open(self.tr_fn, 'r') as fh:
            name = fh.readline()
            if name == "":  # no hits
                return
            name = name.strip()[1:]
            while True:
                line = fh.readline()
                if line == "":
                    break
                if line.startswith("@"):
                    name = line.strip()[1:]
                    continue
                line = line.strip().split(' ')
                data = {x[0]: x[1](y) for x, y in zip(
                    trf_cols, line) if not x[0].startswith("unk")}
                self.trf_lookup[name][data["repeat"]] = data

    def annotate(self, entry, key, new_header):
        """
        Edit the entry if it has a hit
        """
        def edit_entry(repeat, entry, new_header, diff=None):
            # put in the annotations
            entry.translate(new_header)
            entry.info["TRF"] = True
            if diff:
                entry.info["TRFDiff"] = diff
            entry.info["TRFperiod"] = repeat["period"]
            entry.info["TRFcopies"] = repeat["copies"]
            #entry.info["TRFconsize"] = repeat["consize"]
            entry.info["TRFscore"] = repeat["score"]
            entry.info["TRFentropy"] = repeat["entropy"]
            entry.info["TRFrepeat"] = repeat["repeat"]
            return entry

        # 1 - I want trf hits that have srep hits
        if key in self.trf_lookup:
            for repeat in sorted(self.trf_lookup[key].keys(), key=len, reverse=True):
                if repeat in self.srep_lookup[key]:
                    diff = self.trf_lookup[key][repeat]['copies'] - \
                        self.srep_lookup[key][repeat]['copies']
                    repeat = self.trf_lookup[key][repeat]
                    return edit_entry(repeat, entry, new_header, diff)

        # 2 - if there are none, I'll take srep hits
        if key in self.srep_lookup:
            repeat = sorted(
                self.srep_lookup[key].keys(), key=len, reverse=True)[0]
            repeat = self.srep_lookup[key][repeat]
            return edit_entry(repeat, entry, new_header)
        # 3 - otherwise quit
        return entry


def edit_header(header):
    """
    New VCF INFO fields
    """
    header = header.copy()
    header.add_line(('##INFO=<ID=TRF,Number=1,Type=Flag,'
                     'Description="Entry hits a simple repeat region">'))
    header.add_line(('##INFO=<ID=TRFDiff,Number=1,Type=Float,'
                     'Description="Simple repeat copy difference">'))
    header.add_line(('##INFO=<ID=TRFperiod,Number=1,Type=Integer,'
                     'Description="Period size of the repeat">'))
    header.add_line(('##INFO=<ID=TRFcopies,Number=1,Type=Float,'
                     'Description="Number of copies aligned with the consensus pattern">'))
    header.add_line(('##INFO=<ID=TRFscore,Number=1,Type=Integer,'
                     'Description="TRF Alignment score">'))
    header.add_line(('##INFO=<ID=TRFentropy,Number=1,Type=Float,'
                     'Description="TRF Entropy measure">'))
    header.add_line(('##INFO=<ID=TRFrepeat,Number=1,Type=String,'
                     'Description="TRF Repeat found on entry">'))

    return header


def trf_main(cmdargs):
    """ TRF annotation """
    args = parse_args(cmdargs)
    trfshared.args = args
    v = pysam.VariantFile(trfshared.args.input)
    new_header = edit_header(v.header)

    m_regions = truvari.ref_ranges(
        args.reference, chunk_size=int(args.chunk_size * 1e6))
    with multiprocessing.Pool(args.threads, maxtasksperchild=1) as pool:
        chunks = pool.imap_unordered(process_entries, m_regions)
        pool.close()
        with open(args.output, 'w') as fout:
            fout.write(str(new_header))
            for i in chunks:
                fout.write(i[3])
        pool.join()

    logging.info("Finished trf")


if __name__ == '__main__':
    trf_main(sys.argv[1:])
