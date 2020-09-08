"""
Remap VCF'S alleles sequence to the reference to annotate REMAP = novel/tandem/interspersed
Currently only works on Insertions
"""
import sys
import logging
import argparse
import multiprocessing

import pysam
import truvari
from bwapy import BwaAligner

from truvari.annos.grm import cigmatch


class Remap():
    """ Class for remapping annotation """
    def __init__(self, in_vcf, reference, out_vcf="/dev/stdout", min_length=50, threshold=0.8):
        """ The setup """
        self.in_vcf = in_vcf
        self.reference = reference
        self.out_vcf = out_vcf
        self.min_length = min_length
        self.threshold = threshold
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
        self.n_header = header

    def get_end(self, pos, cigar):
        """
        Expand a cigar to get the end position?
        And how much of the query is used?
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
        
        num_hits = 0
        closest_hit = None
        close_dist = None
        for aln in self.aligner.align_seq(seq):
            # Take out the 'same spot' alignment for deletions
            if is_del and aln.rname == entry.chrom and abs(aln.pos - entry.pos) < len(seq):
                continue
            num_hits += 1 
            if aln.rname != entry.chrom:
                continue
            dist = abs(aln.pos - entry.pos)
            if close_dist is None or dist < close_dist:
                end, soft = self.get_end(aln.pos, aln.cigar)
                pct_query = 1 - soft / len(seq)
                if pct_query >= threshold:
                    close_dist = dist
                    closest_hit = aln
                
        if num_hits == 0:
            return "novel"
        elif close_dist and close_dist <= len(seq):
            return "tandem"

        return "interspersed"

    def annotate_entry(self, entry):
        """
        Annotates entries in the vcf and writes to new vcf
        """
        if truvari.entry_size(entry) >= self.min_length:
            entry = truvari.copy_entry(entry, self.n_header)
            entry.info["REMAP"] = self.remap_entry(entry)
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
    def restricted_float(x):
        x = float(x)
        if x < 0.0 or x > 1.0:
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
        return x

    parser = argparse.ArgumentParser(prog="remap", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", default="/dev/stdin",
                        help="Input VCF (%(default)s)")
    parser.add_argument("-r", "--reference", required=True,
                        help="BWA indexed reference")
    parser.add_argument("-o", "--output", default="/dev/stdout",
                        help="Output VCF (%(default)s)")
    parser.add_argument("-m", "--minlength", default=50, type=int,
                        help="Smallest length of allele to remap (%(default)s)")
    parser.add_argument("-t", "--threshold", type=restricted_float, default=.8,
                        help="Threshold for pct of allele covered (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    return args

def remap_main(cmdargs):
    """
    Program's entry point
    """
    args = parse_args(cmdargs)
    anno = Remap(in_vcf=args.input,
                 reference=args.reference,
                 out_vcf=args.output,
                 min_length=args.minlength,
                 threshold=args.threshold)
    anno.annotate_vcf()
    logging.info("Finished Remap Anno")

if __name__ == '__main__':
    remap_main(sys.argv[1:])
