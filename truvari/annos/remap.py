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

from truvari.annos.grm import cigmatch, cig_pctsim, ref_ranges, line_to_entry

#remap_shared = types.SimpleNamespace()

def get_end(pos, cigar):
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


def remap_entry(aligner, entry, threshold=.8):
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
    for aln in aligner.align_seq(seq):
        # Take out the 'same spot' alignment for deletions
        if is_del and aln.rname == entry.chrom and abs(aln.pos - entry.pos) < len(seq):
            continue
        num_hits += 1 
        if aln.rname != entry.chrom:
            continue
        dist = abs(aln.pos - entry.pos)
        if close_dist is None or dist < close_dist:
            end, soft = get_end(aln.pos, aln.cigar)
            pct_query = 1 - soft / len(seq)
            if pct_query >= threshold:
                close_dist = dist
                closest_hit = aln
            
    if num_hits == 0:
        return "novel"
    elif close_dist and close_dist <= len(seq):
        return "tandem"

    return "interspersed"

def edit_header(my_vcf):
    """
    Add INFO for new field to vcf
    """
    header = my_vcf.header.copy()
    header.add_line(('##INFO=<ID=REMAP,Number=1,Type=String,'
                     'Description="Annotation of alt-seq via">'))
    return header

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
    parser.add_argument("-m", "--minlength", default=50, type=int,
                        help="Smallest length of allele to remap (%(default)s)")
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
    aligner = BwaAligner(args.reference, options="-a")
    with pysam.VariantFile(args.input) as fh:
        header = edit_header(fh)
        out = pysam.VariantFile(args.output, 'w', header=header)
        for entry in fh:
            if truvari.entry_size(entry) >= args.minlength:
                entry = truvari.copy_entry(entry, header)
                entry.info["REMAP"] = remap_entry(aligner, entry)
            out.write(entry)
    logging.info("Finished Remap anno")

if __name__ == '__main__':
    remap_main(sys.argv[1:])
