"""
Remap VCF'S alleles sequence to the reference to annotate REMAP

- novel : Allele has no hits in reference
- tandem : Allele's closest hit is within len(allele) bp of the SV's position
- interspersed : Allele's closest hit is not tandem
- partial : Allele only has partial hit(s) less than --threshold

Which alleles and alignments to consider can be altered with:
- --minlength : minimum SV length to considred (50)
- --dist : For deletion SVs, do not consider alignments that hit within Nbp of the SV's position
(a.k.a. alignments back to the source sequence) (10)
- --threshold : Minimum percent of allele's sequence used by alignment to be considered (.8)
"""
import sys
import bisect
import logging
import argparse

import pysam
try:
    from bwapy import BwaAligner
    HASBWALIB = True
except (OSError,  ModuleNotFoundError):
    HASBWALIB = False

import truvari
from truvari.annotations.grm import cigmatch


class Remap():
    """ Class for remapping annotation """

    def __init__(self, in_vcf, reference, out_vcf="/dev/stdout", min_length=50,
                 threshold=0.8, min_distance=10, anno_hits=0):
        """ The setup """
        self.in_vcf = in_vcf
        self.reference = reference
        self.out_vcf = out_vcf
        self.min_length = min_length
        self.threshold = threshold
        self.min_distance = min_distance
        self.anno_hits = anno_hits
        self.n_header = None
        self.aligner = BwaAligner(self.reference, options="-a")

    def edit_header(self, header=None):
        """
        Edits and holds on to the header
        """
        if header is None:
            with pysam.VariantFile(self.in_vcf, 'r') as fh:
                header = fh.header.copy()
        header.add_line(('##INFO=<ID=REMAP,Number=1,Type=String,'
                         'Description="Annotation of alt-seq remapping">'))
        if self.anno_hits:
            header.add_line(('##INFO=<ID=REMAPHits,Number=.,Type=String,'
                             'Description="List of chr:start-end of hits">'))

        self.n_header = header

    def get_end(self, pos, cigar):  # pylint: disable=no-self-use
        """
        Expand a cigar to get the end position and how many bases are clipped
        """
        soft_bases = 0
        for i in cigmatch.findall(cigar):
            if i[-1] == "S":
                soft_bases += int(i[:-1])
            elif i[-1] in ["M", "D"]:
                pos += int(i[:-1])
        return pos, soft_bases

    def remap_entry(self, entry, threshold=.8):
        """
        Map a sequence and return the information from it
        """
        is_del = truvari.entry_variant_type(entry) == "DEL"
        if is_del:
            seq = str(entry.ref)
        else:
            seq = entry.alts[0]
        all_hits = []
        num_hits = 0
        partial_hits = 0
        close_dist = None
        for aln in self.aligner.align_seq(seq):
            # Take out the 'same spot' alignment for deletions
            dist = abs(aln.pos - entry.pos)
            if is_del and aln.rname == entry.chrom and dist < self.min_distance:
                continue

            # Filter hits below threshold
            end, soft = self.get_end(aln.pos, aln.cigar)
            seq_len = len(seq)
            pct_query = (seq_len - soft) / seq_len
            if pct_query < threshold:
                partial_hits += 1
                continue
            hit = f"{aln.rname}:{aln.pos}-{end}.{int(pct_query*100)}"
            bisect.insort(all_hits, (pct_query, hit))
            num_hits += 1
            if aln.rname != entry.chrom:
                continue
            if close_dist is None or dist < close_dist:
                close_dist = dist

        if num_hits == 0 and partial_hits == 0:
            return "novel", all_hits
        if close_dist and close_dist <= len(seq):
            return "tandem", all_hits
        if num_hits == 0 and partial_hits != 0:
            return "partial", all_hits

        return "interspersed", all_hits

    def annotate_entry(self, entry):
        """
        Annotates entries in the vcf and writes to new vcf
        """
        if truvari.entry_size(entry) >= self.min_length:
            entry.translate(self.n_header)
            remap, hits = self.remap_entry(entry)
            entry.info["REMAP"] = remap
            if self.anno_hits and hits:
                entry.info["REMAPHits"] = [_[1] for _ in hits[-self.anno_hits:]]
        return entry

    def annotate_vcf(self):
        """
        Annotates the vcf and writes to new vcf
        """
        with pysam.VariantFile(self.in_vcf) as fh:
            self.edit_header(fh.header.copy())
            out = pysam.VariantFile(self.out_vcf, 'w', header=self.n_header)
            for entry in fh:
                entry = self.annotate_entry(entry)
                out.write(entry)


def parse_args(args):
    """
    Argument parsing
    """
    parser = argparse.ArgumentParser(prog="remap", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", default="/dev/stdin",
                        help="Input VCF (%(default)s)")
    parser.add_argument("-r", "--reference", required=True,
                        help="BWA indexed reference")
    parser.add_argument("-o", "--output", default="/dev/stdout",
                        help="Output VCF (%(default)s)")
    parser.add_argument("-m", "--minlength", default=50, type=truvari.restricted_int,
                        help="Smallest length of allele to remap (%(default)s)")
    parser.add_argument("-t", "--threshold", type=truvari.restricted_float, default=.8,
                        help="Threshold for pct of allele covered to consider hit (%(default)s)")
    parser.add_argument("-d", "--dist", type=truvari.restricted_int, default=10,
                        help=("Minimum distance an alignment must be from a DEL's "
                              "position to be considered (%(default)s))"))
    parser.add_argument("-H", "--hits", type=truvari.restricted_int, default=0,
                        help="Report top hits as chr:start-end.pct (max %(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    return args


def remap_main(cmdargs):
    """
    Program's entry point
    """
    if not HASBWALIB:
        logging.error("bwapy isn't available on this machine")
        sys.exit(1)
    args = parse_args(cmdargs)
    anno = Remap(in_vcf=args.input,
                 reference=args.reference,
                 out_vcf=args.output,
                 min_length=args.minlength,
                 threshold=args.threshold,
                 anno_hits=args.hits)
    anno.annotate_vcf()
    logging.info("Finished remap")
