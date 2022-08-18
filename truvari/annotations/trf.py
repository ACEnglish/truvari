"""
Intersect vcf with reference simple repeats and report
how the an alternate allele affects the copies using TRF
"""
import os
import sys
import json
import types
import shutil
import logging
import argparse
import multiprocessing
from io import StringIO
from functools import cmp_to_key

import pysam
import truvari

trfshared = types.SimpleNamespace()

try:
    from setproctitle import setproctitle  # pylint: disable=import-error,useless-suppression
except ModuleNotFoundError:
    def setproctitle(_):
        """ dummy function """
        return

def compare_scores(a, b):
    """
    sort annotations
    """
    # most amount of SV covered
    ret = 0
    if a['ovl_pct'] > b['ovl_pct']:
        ret = 1
    elif a['ovl_pct'] < b['ovl_pct']:
        ret = -1
    elif a['score'] > b['score']:
        ret = 1
    elif a['score'] < b['score']:
        ret = -1
    else:
        aspan = a['end'] - a['start']
        bspan = b['end'] - b['start']
        if aspan > bspan:
            ret = 1
        elif aspan < bspan:
            ret = -1
    return ret
score_sorter = cmp_to_key(compare_scores)

class TRFAnno():
    """
    Class for trf annotation
    Operates on a single TRF region across multiple TRF annotations
    """

    def __init__(self, region, reference, motif_similarity=0.90):
        """ setup """
        self.region = region
        self.reference = reference
        self.motif_similarity = motif_similarity
        self.known_motifs = {_['repeat']:_['copies'] for _ in self.region['annos']}

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

    def get_anno_span(self, entry):
        """
        Get the subsequence positions of self.reference that maximally cover the entry
        For example, an INS inside this region may only intersect some of the self.region['annos'],
        Iterate those and get the min_start/max_end of those annos
        If the entry doesn't hit any annnotations in the region, return 0, len(self.reference)
        """
        min_start = self.region['end']
        max_end = self.region['start']
        for i in self.region['annos']:
            if i['start'] <= entry.start < i['end']:
                min_start = min(min_start, i['start'])
            if i['start'] <= entry.stop < i['end']:
                max_end = max(max_end, i['end'])
        # variant's start is not inside an annotation. set it to the region start
        if min_start == self.region['end']:
            min_start = self.region['start']
        if max_end == self.region['start']:
            max_end = self.region['end']
        return min_start, max_end

    def make_seq(self, entry, svtype, substr=None):
        """
        Make the haplotype sequence
        use r_start/r_end to make sequence a subsequence
        """
        # variant position relative to this region
        r_start = entry.start - self.region['start']
        r_end = entry.stop - self.region['start']
        seq = self.reference
        r_start = entry.start - self.region['start']
        r_end = entry.stop - self.region['start']
        if substr:
            seq = seq[substr[0]:substr[1]]
            r_start += substr[0]
            r_end += substr[0]
        up_seq = seq[:r_start]
        dn_seq = seq[r_end:]
        if svtype == "INS":
            m_seq = up_seq + entry.alts[0] + dn_seq
        elif svtype == "DEL":
            m_seq = up_seq + dn_seq
        else:
            logging.critical("Can only consider entries with 'SVTYPE' INS/DEL")
            sys.exit(1)
        return m_seq

    def score_annotation(self, var_start, var_end, anno):
        """
        Scores the annotation. Addes fields in place.
        if is_new, we calculate the diff
        """
        ovl_pct = truvari.overlap_percent(var_start, var_end, anno['start'], anno['end'])
        # has to have overlap
        if ovl_pct <= 0:
            return None
        anno['ovl_pct'] = ovl_pct

        if anno['repeat'] in self.known_motifs:
            anno['diff'] = anno['copies'] - self.known_motifs[anno['repeat']]
        else:
            best_score = 0
            best_pos = None
            motif_len = len(anno['repeat'])
            for known in self.region['annos']:
                if truvari.sizesim(len(known['repeat']), motif_len)[0] < self.motif_similarity:
                    continue
                sq = truvari.unroll_compare(known['repeat'], anno['repeat'], anno['start'] - known['start'])
                if sq >= self.motif_similarity and sq > best_score:
                    best_score = sq
                    best_pos = known
   
            if best_pos:
                anno['diff'] = anno['copies'] - best_pos['copies']
                anno['orig_repeat'] = anno['repeat']
                anno['repeat'] = best_pos['repeat']
            else:
                anno['diff'] = 0

        return anno

    def del_annotate(self, entry, svlen, score_filter=True):
        """
        Annotate a deletion
        """
        scores = []
        for anno in self.region['annos']:
            ovl_pct = truvari.overlap_percent(entry.start, entry.stop, anno['start'], anno['end'])
            if ovl_pct == 0:
                continue
            m_sc = dict(anno)
            m_sc['ovl_pct'] = ovl_pct
            m_sc['diff'] = - (ovl_pct * svlen) / anno['period']
            scores.append(m_sc)

        scores.sort(reverse=True, key=score_sorter)
        if score_filter and scores:
            return scores[0]
        return scores

    def ins_annotate(self, entry, annos=None, score_filter=True):
        """
        Annotate an insertion
        If candidate insertions are made elsewhere (e.g. from batching),
        use those. otherwise this method will run_trf itself
        """
        if annos is None:
            seq = self.make_seq(entry, 'INS')
            fa_fn = truvari.make_temp_filename(suffix='.fa')
            with open(fa_fn, 'w') as fout:
                fout.write(f">key\n{seq}\n")
            annos = run_trf(fa_fn)
            if annos:
                annos = annos['key']
            # Shift back to genomic coordinates
            for i in annos:
                i['start'] += self.region['start']
                i['end'] += self.region['end']
                

        scores = []
        for anno in annos:
            m_sc = self.score_annotation(entry.start, entry.stop, anno)
            if m_sc:
                scores.append(m_sc)
        scores.sort(reverse=True, key=score_sorter)
        if score_filter and scores:
            return scores[0]
        return scores

    def annotate(self, entry, score_filter=True):
        """
        Figure out the hit and return
        """
        svtype = truvari.entry_variant_type(entry)
        sz = truvari.entry_size(entry)
        repeat = []
        if svtype == 'DEL':
            repeat = self.del_annotate(entry, sz, score_filter)
        elif svtype == 'INS':
            repeat = self.ins_annotate(entry, score_filter=score_filter)
        return repeat

    def translate_coords(self, annos):
        """
        Given a trf result, translate the coordantes back into this region
        """
        for anno in annos:
            anno['start'] += self.region['start']
            anno['end'] += self.region['start']

def parse_trf_output(fn):
    """
    Parse the outputs from TRF
    Returns a list of hits
    """
    # The columns I don't want are set to type None
    trf_cols = [("start", int),
                ("end", int),
                ("period", int),
                ("copies", float),
                ("consize", int),
                ("pctmat", None), # int
                ("pctindel", None), # int
                ("score", int),
                ("A", None), # int
                ("C", None), # int
                ("G", None), # int
                ("T",  None), # int
                ("entropy", float),
                ("repeat", str),
                ("unk1", None), # str
                ("unk2", None), # str
                ("unk3", None) # str
            ]

    ret = {}
    with open(fn, 'r') as fh:
        var_key = None
        for line in fh:
            line = line.strip()
            if line.startswith("@"):
                var_key = line[1:]
                ret[var_key] = []
                continue

            data = {x[0]: x[1](y) for x, y in zip(trf_cols, line.split(' ')) if x[1]}
            # correction to 0-based
            data['start'] -= 1
            ret[var_key].append(data)

    return ret

def run_trf(fa_fn, executable="trf409.linux64",
            trf_params="2 7 7 80 10 50 500 -m -f -h -d -ngs"):
    """
    Given a fasta file, run TRF and return result
    """
    if "-ngs" not in trf_params:
        trf_params = trf_params + " -ngs "

    tr_fn = fa_fn + '.txt'
    cmd = f"{executable} {fa_fn} {trf_params} > {tr_fn}"
    ret = truvari.cmd_exe(cmd)
    if ret.ret_code != 0:
        logging.error("Couldn't run trf. Check Parameters")
        logging.error(cmd)
        logging.error(str(ret))
        return {}

    return parse_trf_output(tr_fn)


def process_ref_region(region):
    """
    Process a section of the reference.
    Tries to run TRF only once
    """
    logging.debug(f"Starting region {region['chrom']}:{region['start']}-{region['end']}")
    setproctitle(f"trf {region['chrom']}:{region['start']}-{region['end']}")

    vcf = pysam.VariantFile(trfshared.args.input)
    to_consider = []
    try:
        for entry in vcf.fetch(region['chrom'], region['start'], region['end']):
            # Variants must start within region
            if not (entry.start >= region['start'] and entry.start < region['end']):
                continue
            to_consider.append(entry)
    except ValueError as e:
        logging.debug("Skipping VCF fetch %s", e)

    """
    1 - I'll need to pull all the annotations . So we can use tabix to just pull them out for this ref region
    Every Region will need its own TRFAnno
    anno_regions = [] (list of TRFAnnos)
    I'll build these as we traverse the regions

    If it's a DEL, we can process it immediately and write

    If it's an INS, we can build the sequence (will need keys.. which point to entry and TRFAnnos)
        (I also want to 'shorten' the sequences we're looking at to as few bases as possible - 
        So Don't annotate the entire reference sequence, trf on just the annotation sub-region

    Then we can submit all the INS at once
    Then we can...
        run_batch_trf
        Which will create from entries.. we have options here
        
        This will return the dictionary from parse_trf_output
        And this will make keys of the position in the list
    Then we can iterate the annotations per position on the list
    for key, (entry, TRFAnno) in enumerate(all_ins_were_looking_at):
        annos = ret from run_batch_trf[key]
        TRFAnno.annotate(entry, annos)
        write
    Done
    
    So.. this will reduce the overhead of creating lots of processes.
    And if I can get the coordinate shifting to so that INS are only over the max span, that'd 
    help things out as well

    TODO:
    A) need to try an shorten the coordinates to only the maximum intersecting span
        put into the key a way to identify the variant (int position will work) and the absolute start of the subseq
    I) Add a --region option (either a string of a file) so that we don't always have to iterate the whole genome
    J) I could add (just do it per-region and extract the variants) then process_tr_region doesn't exactly go to waste
    COMPLETED:
    B) Extract run_trf from TRFAnno (It shouldn't actually need anything inside TRFAnno
    C) Make the ins_annotation thing write its own fasta and then call the extracted run_trf
    D) The TRFAnno.annotate shold return the annotation (repeat annotation dict)
    E) Add TRFAnno.edit_entry(entry, annotation) - static method, probably
    F) main doesn't need to do all this crazy pre-reading of the VCF stuff. each chunk will work
    G) For chunking only the START needs to be in the chunk-region
    H) Need to get back the reference chunking stuff

    """

def process_tr_region(region):
    """
    Process vcf lines from a tr reference section
    """
    logging.debug(f"Starting region {region['chrom']}:{region['start']}-{region['end']}")
    setproctitle(f"trf {region['chrom']}:{region['start']}-{region['end']}")

    ref = pysam.FastaFile(trfshared.args.reference)
    ref_seq = ref.fetch(region['chrom'], region['start'], region['end'])
    tanno = TRFAnno(region, ref_seq, trfshared.args.motif_similarity)
    vcf = pysam.VariantFile(trfshared.args.input)
    new_header = edit_header(vcf.header)
    out = StringIO()

    fa_fn = truvari.make_temp_filename(suffix='.fa')
    batch = []
    batch_size = 0
    try:
        m_fetch = vcf.fetch(region['chrom'], region['start'], region['end'])
    except ValueError as e:
        logging.debug("Skipping VCF fetch %s", e)
    
    with open(fa_fn, 'w') as fa_out:
        for entry in m_fetch:
            # Variants must be entirely contained within region
            if not (entry.start >= region['start'] and entry.stop < region['end']):
                continue
            entry.translate(new_header)
            svtype = truvari.entry_variant_type(entry)
            svlen = truvari.entry_size(entry)
            if svlen < trfshared.args.min_length or svtype not in ["DEL", "INS"]:
                edit_entry(entry, None)
                out.write(str(entry))
                continue
    
            if svtype == 'DEL':
                m_anno = tanno.del_annotate(entry, svlen)
                edit_entry(entry, m_anno)
                out.write(str(entry))
            elif svtype == 'INS':
                # translate coords to subsequence of region
                # This is broken. And it makes sense we get more consistent
                # Annotations from the same reference span
                # r_start, r_end = tanno.get_anno_span(entry)
                # sub_start = r_start - region['start']
                # sub_end = r_end - region['start']
                seq = tanno.make_seq(entry, 'INS')
                fa_out.write(f">{batch_size}\n{seq}\n")
                #batch.append((entry, r_start))
                batch.append(entry)
                batch_size += 1
    
    if batch:
        annotations = run_trf(fa_fn, trfshared.args.executable, trfshared.args.trf_params)
        for idx, entry in enumerate(batch):
            key = str(idx)
            m_anno = None
            if key in annotations:
                # translate coords back
                for i in annotations[key]:
                    i['start'] += region['start']
                    i['end'] += region['start']
                m_anno = tanno.ins_annotate(entry, annotations[key])
            edit_entry(entry, m_anno)
            out.write(str(entry))

    out.seek(0)
    setproctitle(f"trf {region['chrom']}:{region['start']}-{region['end']}")
    logging.debug(f"Done region {region['chrom']}:{region['start']}-{region['end']}")
    return out.read()

def iter_tr_regions(fn, coords=None):
    """
    Read a repeats file with structure chrom, start, end, annotations.json
    returns generator of dicts
    """
    for line in truvari.opt_gz_open(fn):
        chrom, start, end, annos = line.strip().split('\t')
        start = int(start)
        end = int(end)
        annos = json.loads(annos)
        yield {'chrom': chrom,
               'start': start,
               'end': end,
               'annos': annos}

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

def edit_entry(entry, repeat):
    """
    places the INFO fields into the entry
    assumes entry has been translated to a vcf header holding the correct fields
    """
    entry.info['TRF'] = True
    if repeat:
        entry.info["TRFovl"] = round(repeat['ovl_pct'], 3)
        entry.info["TRFdiff"] = round(repeat['diff'], 1)
        entry.info["TRFperiod"] = repeat["period"]
        entry.info["TRFcopies"] = repeat["copies"]
        entry.info["TRFscore"] = repeat["score"]
        entry.info["TRFentropy"] = repeat["entropy"]
        entry.info["TRFrepeat"] = repeat["repeat"]

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
    parser.add_argument("-s", "--motif-similarity", type=truvari.restricted_float, default=0.90,
                        help="Motif similarity score (%(default)s)")
    parser.add_argument("-m", "--min-length", type=truvari.restricted_int, default=50,
                        help="Minimum size of entry to annotate (%(default)s)")
    parser.add_argument("-M", "--max-length", type=truvari.restricted_int, default=10000,
                        help="Maximum size of sequence to run through trf (%(default)s)")
    parser.add_argument("-R", "--regions-only", action='store_true',
                        help="Only write variants within --repeats regions (%(default)s)")
    parser.add_argument("-C", "--chunk-size", type=int, default=1,
                            help="Size (in mbs) of reference chunks for parallelization (%(default)s)")
    parser.add_argument("-t", "--threads", type=truvari.restricted_int, default=multiprocessing.cpu_count(),
                        help="Number of threads to use (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")

    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    return args

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

def trf_single_main(cmdargs):
    """ TRF annotation but single threaded (for debugging/profiling)"""
    args = parse_args(cmdargs)
    check_params(args)
    trfshared.args = args

    m_regions = iter_tr_regions(args.repeats)

    vcf = pysam.VariantFile(trfshared.args.input)
    new_header = edit_header(vcf.header)

    with open(args.output, 'w') as fout:
        fout.write(str(new_header))
        for i in m_regions:
            fout.write(process_tr_region(i))
    logging.info("Finished trf")

def trf_main(cmdargs):
    """ TRF annotation """
    args = parse_args(cmdargs)
    check_params(args)
    trfshared.args = args
    
    m_regions = None
    m_process = None
    if args.regions_only:
        m_regions = iter_tr_regions(args.repeats)
        m_process = process_tr_region
    else:
        m_regions = truvari.ref_ranges(args.reference, chunk_size=int(args.chunk_size * 1e6))
        m_process = process_ref_region
    

    vcf = pysam.VariantFile(trfshared.args.input)
    new_header = edit_header(vcf.header)

    with multiprocessing.Pool(args.threads, maxtasksperchild=1) as pool:
        chunks = pool.imap_unordered(m_process, m_regions)
        pool.close()
        with open(args.output, 'w') as fout:
            fout.write(str(new_header))
            for i in chunks:
                fout.write(i)
        pool.join()

    logging.info("Finished trf")


if __name__ == '__main__':
    trf_main(sys.argv[1:])
