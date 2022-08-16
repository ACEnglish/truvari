"""
Intersect vcf with reference simple repeats and report
how the an alternate allele affects the copies using TRF
"""
import os
import sys
import json
import types
import shutil
import decimal
import logging
import argparse
import tempfile
import multiprocessing
from io import StringIO
from functools import cmp_to_key
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

def parse_trf_output(fn):
    """
    Parse the outputs from TRF
    Returns a list of hits
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
                ("unk1", str),
                ("unk2", str),
                ("unk3", str)]
    ret = defaultdict(list)
    with open(fn, 'r') as fh:
        var_key = None
        while True:
            line = fh.readline().strip()
            if line == "":
                break
            if line.startswith("@"):
                var_key = line[1:]
                continue
            data = {x[0]: x[1](y) for x, y in zip(trf_cols, line.split(' '))}
            # 0-based correction
            data['start'] -= 1
            ret[var_key].append(data)
    return dict(ret)

def compare_scores(a, b):
    """
    comparator for scoring
    """
    # Prefer reference first
    if 'diff' in a[2] and 'diff' not in b[2]:
        return 1
    if 'diff' not in a[2] and 'diff' in b[2]:
        return -1
    # Overlap
    if a[0] > b[0]:
        return 1
    if a[0] < b[0]:
        return -1
    # Score
    if a[1] > b[1]:
        return 1
    if a[1] < b[1]:
        return -1
    # Take the longer motif
    if len(a[2]["repeat"]) > len(b[2]["repeat"]):
        return 1
    if len(a[2]["repeat"]) < len(b[2]["repeat"]):
        return -1
    return 0
score_sorter = cmp_to_key(compare_scores)

class TRFAnno():
    """
    Class for trf annotation
    Operates on a single TRF region across multiple TRF annotations
    """

    def __init__(self, chrom, start, end, reference, repeats,
                 min_length=50, max_length=10000,
                 executable="trf409.linux64",
                 trf_params="2 7 7 80 10 50 500 -m -f -h -d -ngs",
                 tmpdir=None):
        """ setup """
        self.region_chrom = chrom
        self.region_start = start
        self.region_end = end
        self.reference = reference
        self.repeats = repeats
        self.min_length = min_length
        self.max_length = max_length
        self.executable = executable
        if "-ngs" not in trf_params:
            trf_params = trf_params + " -ngs "
        self.trf_params = trf_params

        # Where we write the fasta entries
        self.fa_fn = truvari.make_temp_filename(tmpdir, '.fa')
        self.tr_fn = self.fa_fn + '.txt'
        # Populated later
        self.unfilt_annotations = None # 1-to-N entrykey to trf results
        self.annotations = None # 1-to-1 entrykey to trf result
        self.known_motifs = {_['repeat']:_['copies'] for _ in self.repeats}

    def entry_to_key(self, entry):
        """
        VCF entries to names for the fa and header lines in the tr output
        returns the key and the entry's size
        """
        sz = truvari.entry_size(entry)
        svtype = truvari.entry_variant_type(entry)
        o_sz = sz if svtype == 'INS' else 0 # span of variant in alt-seq
        key = f"{entry.chrom}:{entry.start}:{entry.stop}:{o_sz}:{hash(entry.ref)}:{hash(entry.alts[0])}"
        return key, sz, svtype

    def make_seq(self, entry, svtype):
        """
        Make the haplotype sequence
        """
        if svtype == "INS":
            m_seq = self.reference[:entry.start - self.region_start] + \
                entry.alts[0] + self.reference[entry.stop - self.region_start:]
        elif svtype == "DEL":
            m_start = max(self.region_start, entry.start)
            m_end = min(self.region_end, entry.stop)
            m_seq = self.reference[:m_start - self.region_start] + self.reference[m_end - self.region_start:]
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
                key, sz, svtype = self.entry_to_key(entry)
                if sz < self.min_length:
                    continue
                seq = self.make_seq(entry, svtype)
                if self.min_length <= len(seq) <= self.max_length:
                    n_seqs += 1
                    fout.write(f">{key}\n{seq}\n")
        if not n_seqs:
            self.annotations = {}
            return

        ret = truvari.cmd_exe(
            f"{self.executable} {self.fa_fn} {self.trf_params} > {self.tr_fn}")
        if ret.ret_code != 0:
            logging.error("Couldn't run trf. Check Parameters")
            logging.error(f"{self.executable} {self.fa_fn} {self.trf_params} > {self.tr_fn}")
            logging.error(str(ret))
            self.annotations = {}
            return

        self.unfilt_annotations = parse_trf_output(self.tr_fn)
        self.filter_annotations()

    def filter_annotations(self):
        """
        Pick the best annotation for every entry in self.unfilt_annotations
        places results into self.annotations
        I've put this into its own method because eventually, maybe, we can have options on how to choose
        """
        self.annotations = {}
        for var_key in self.unfilt_annotations:
            var_start, var_end, var_len = var_key.split(":")[1:4]
            var_start = int(var_start)
            var_end = int(var_end)
            var_len = int(var_len)
            scores = []
            for m_anno in self.unfilt_annotations[var_key]:
                m_start = var_start - self.region_start
                m_end = m_start + var_len
                m_sc = self.score_annotation(m_start, m_end, m_anno, True)
                if m_sc:
                    logging.debug("made trf score %s", str(m_anno))
                    scores.append(m_sc)
            scores.sort(reverse=True, key=score_sorter)
            if scores:
                self.annotations[var_key] = scores[0][-1]

    def score_annotation(self, var_start, var_end, anno, is_new=False):
        """
        Scores the annotation. Addes fields in place.
        if is_new, we calculate the diff
        """
        ovl_pct = truvari.overlap_percent(var_start, var_end, anno['start'], anno['end'])
        logging.debug('ovl')
        logging.debug(ovl_pct)
        logging.debug("%d %d - %d %d", var_start, var_end, anno['start'], anno['end'])
        # has to hit
        if ovl_pct <= 0:
            return None
        if is_new and anno['repeat'] in self.known_motifs:
            logging.debug('in motifs? %s', anno['repeat'])
            anno['diff'] = anno['copies'] - self.known_motifs[anno['repeat']]
        anno['ovl_pct'] = ovl_pct
        return ovl_pct, anno['score'], anno

    def annotate(self, entry):
        """
        Edit the entry if it has a hit
        """
        def edit_entry(entry, repeat):
            if 'ovl_pct' in repeat:
                entry.info["TRFovl"] = round(repeat['ovl_pct'], 3)
            if 'diff' in repeat:
                entry.info["TRFdiff"] = repeat['diff']
            entry.info["TRFperiod"] = repeat["period"]
            entry.info["TRFcopies"] = repeat["copies"]
            entry.info["TRFscore"] = repeat["score"]
            entry.info["TRFentropy"] = repeat["entropy"]
            entry.info["TRFrepeat"] = repeat["repeat"]

        entry.info["TRF"] = True
        key, sz, _ = self.entry_to_key(entry)
        if sz < self.min_length:
            return
        if key in self.annotations:
            logging.debug("okay, we're inside")
            edit_entry(entry, self.annotations[key])
            return

        # pick best reference annotation
        scores = []
        for m_anno in self.repeats:
            m_sc = self.score_annotation(entry.start, entry.stop, dict(m_anno))
            if m_sc:
                logging.debug("made ref score")
                scores.append(m_sc)
        scores.sort(reverse=True, key=score_sorter)
        if scores:
            edit_entry(entry, scores[0][-1])

    def cleanup(self):
        """
        remove temporary files
        """
        try:
            os.remove(self.fa_fn)
        except FileNotFoundError:
            pass
        try:
            os.remove(self.tr_fn)
        except FileNotFoundError:
            pass

def process_tr_region(region):
    """
    Process vcf lines from a tr reference section
    """
    r_chrom, r_start, r_end, r_annos = region
    logging.debug(f"Starting region {r_chrom}:{r_start}-{r_end}")
    setproctitle(f"trf {r_chrom}:{r_start}-{r_end}")
    vcf = pysam.VariantFile(trfshared.args.input)
    to_consider = []
    try:
        for entry in vcf.fetch(r_chrom, r_start, r_end):
            # Variants must be entirely contained within region
            if not (entry.start >= r_start and entry.stop < r_end):
                continue
            to_consider.append(entry)
    except ValueError as e:
        logging.debug("Skipping VCF fetch %s", e)

    # no variants, so nothing to do
    if not to_consider:
        return ""

    ref_seq = pysam.FastaFile(trfshared.args.reference).fetch(r_chrom, r_start, r_end)
    tanno = TRFAnno(r_chrom, r_start, r_end, ref_seq, r_annos,
                    min_length = trfshared.args.min_length,
                    max_length = trfshared.args.max_length,
                    executable=trfshared.args.executable,
                    trf_params=trfshared.args.trf_params,
                    tmpdir=trfshared.args.tmpdir)
    tanno.run_trf(to_consider)

    new_header = edit_header(vcf.header)
    out = StringIO()
    for entry in to_consider:
        entry.translate(new_header)
        tanno.annotate(entry)
        out.write(str(entry))
    out.seek(0)
    # cleanup after yourself
    tanno.cleanup()
    setproctitle(f"trf done {r_chrom}:{r_start}-{r_end}")
    logging.debug(f"Done region {r_chrom}:{r_start}-{r_end}")
    return out.read()


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
    parser.add_argument("-r", "--repeats", type=str, required=True,
                        help="Reference repeat annotations")
    parser.add_argument("-f", "--reference", type=str, required=True,
                        help="Reference fasta file")
    parser.add_argument("-m", "--min-length", type=truvari.restricted_int, default=50,
                        help="Minimum size of entry to annotate (%(default)s)")
    parser.add_argument("-M", "--max-length", type=truvari.restricted_int, default=10000,
                        help="Maximum size of sequence to run through trf (%(default)s)")
    parser.add_argument("-t", "--threads", type=truvari.restricted_int, default=multiprocessing.cpu_count(),
                        help="Number of threads to use (%(default)s)")
    parser.add_argument("--tmpdir", type=str, default=None,
                        help="Temporary directory")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")

    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    if args.tmpdir and not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
    return args



def edit_header(header):
    """
    New VCF INFO fields
    """
    header = header.copy()
    header.add_line(('##INFO=<ID=TRF,Number=1,Type=Flag,'
                     'Description="Entry hits a simple repeat region">'))
    header.add_line(('##INFO=<ID=TRFdiff,Number=1,Type=Float,'
                     'Description="ALT TR copy difference from reference">'))
    header.add_line(('##INFO=<ID=TRFperiod,Number=1,Type=Integer,'
                     'Description="Period size of the repeat">'))
    header.add_line(('##INFO=<ID=TRFcopies,Number=1,Type=Float,'
                     'Description="Number of copies aligned with the consensus pattern">'))
    header.add_line(('##INFO=<ID=TRFscore,Number=1,Type=Integer,'
                     'Description="Alignment score">'))
    header.add_line(('##INFO=<ID=TRFentropy,Number=1,Type=Float,'
                     'Description="Entropy measure">'))
    header.add_line(('##INFO=<ID=TRFrepeat,Number=1,Type=String,'
                     'Description="Repeat motif">'))
    header.add_line(('##INFO=<ID=TRFovl,Number=1,Type=Float,'
                     'Description="Percent of ALT covered by TRF annotation">'))
    return header


def check_params(args):
    """
    Ensure the files are compressed/indexed
    """
    check_fail = False
    if not os.path.exists(args.input):
        logging.error(f"{args.input} doesn't exit")
        check_fail = True
    if not args.input.endswith((".vcf.gz", ".bcf.gz")):
        logging.error(f"{args.input} isn't compressed vcf")
        check_fail = True
    if not os.path.exists(args.input + '.tbi') and not os.path.exists(args.input + '.csi'):
        logging.error(f"{args.input}[.tbi|.csi] doesn't exit")
        check_fail = True
    if not args.repeats.endswith(".bed.gz"):
        logging.error(f"{args.repeats} isn't compressed bed")
        check_fail = True
    if not os.path.exists(args.repeats + '.tbi'):
        logging.error(f"{args.repeats}.tbi doesn't exit")
        check_fail = True
    if not shutil.which(args.executable):
        logging.error(f"{args.executable} not found in path")
        check_fail = True
    if check_fail:
        logging.error("Please fix parameters")
        sys.exit(1)

def iter_tr_regions(fn):
    """
    Read a repeats file with structure chrom, start, end, annotations.json
    returns generator
    """
    for line in truvari.opt_gz_open(fn):
        chrom, start, end, annos = line.strip().split('\t')
        start = int(start)
        end = int(end)
        annos = json.loads(annos)
        yield chrom, start, end, annos

def trf_main(cmdargs):
    """ TRF annotation """
    args = parse_args(cmdargs)
    check_params(args)
    trfshared.args = args
    vcf = pysam.VariantFile(trfshared.args.input)
    new_header = edit_header(vcf.header)

    m_lookup, _ = truvari.build_anno_tree(args.repeats)
    m_regions = iter_tr_regions(args.repeats)

    with multiprocessing.Pool(args.threads, maxtasksperchild=1) as pool:
        chunks = pool.imap_unordered(process_tr_region, m_regions)
        pool.close()
        with open(args.output, 'w') as fout:
            fout.write(str(new_header))
            # Write variants not considered
            for entry in vcf:
                hits = m_lookup[entry.chrom][entry.start:entry.stop]
                has_hit = False
                for i in hits:
                    if entry.start >= i.begin and entry.stop < i.end:
                        has_hit = True
                        break
                if not has_hit:
                    fout.write(str(entry))
            # Now collect the others
            for i in chunks:
                fout.write(i)
        pool.join()

    logging.info("Finished trf")


if __name__ == '__main__':
    trf_main(sys.argv[1:])
