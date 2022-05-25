"""
Segmentation: Normalization of SVs into disjointed genomic regions
"""
import logging
import argparse
from collections import defaultdict

import pysam
import numpy as np
from intervaltree import IntervalTree

import truvari


def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="segment", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("vcf", metavar="IN",
                        help="VCF to parse")
    parser.add_argument("output", metavar="OUT",
                        help="Output VCF")
    # parser.add_argument("-m", "--min", default=10, type=int,
    # help="Minimum span of variants to segment")
    # parser.add_argument("--alter", action="store_true",
    # help="Add SEG Format field to all variants (false)")
    args = parser.parse_args(args)
    truvari.setup_logging()
    return args


def edit_header(vcf):
    """
    Add new vcf header lines
    """
    header = vcf.header.copy()
    header.add_line(('##INFO=<ID=SEGCNT,Number=1,Type=Integer,'
                     'Definition="Number of SVs overlapping this segment">'))
    header.add_line(('##FORMAT=<ID=SEG,Number=1,Type=Integer,'
                     'Definition="Segmentation\'s disjoint region coverage (may be > ploidy)">'))
    return header


def segment_main(args):
    """
    Main entry point for running Segmentation
    """
    args = parse_args(args)

    vcf = pysam.VariantFile(args.vcf)
    header = edit_header(vcf)
    out = pysam.VariantFile(args.output, 'w', header=header)
    tree = defaultdict(IntervalTree)
    # Assume 0 for anything not HET/HOM
    gtcnt = defaultdict(int)
    gtcnt["HET"] = 1
    gtcnt["HOM"] = 2

    # needs to be split by chrom
    for entry in vcf:
        if "SVTYPE" not in entry.info or entry.info["SVTYPE"] != "DEL":
            out.write(entry)
            continue
        data = [gtcnt[truvari.get_gt(x["GT"]).name]
                for x in entry.samples.values()]
        tree[entry.chrom].addi(entry.start, entry.stop, data=([entry], data))

    def reducer(x, y):
        return (x[0] + y[0], np.array(x[1], dtype=np.int8) + np.array(y[1], dtype=np.int8))
    # This should probably be its own method so that as I'm splitting by chrom, I don't have to
    # hold everything in memory
    for chrom, chr_tree in tree.items():
        logging.info(f"Segmenting {chrom}")
        chr_tree.split_overlaps()
        chr_tree.merge_overlaps(data_reducer=reducer)

        for reg in sorted(chr_tree):
            segcnt = len(reg.data[0])
            entry = reg.data[0][0]
            # Unaltered
            if segcnt == 1 and entry.start == reg.begin and entry.stop == reg.end:
                out.write(entry)
                continue
            new_entry = out.new_record()
            new_entry.chrom = entry.chrom
            new_entry.start = reg.begin
            new_entry.stop = reg.end
            # assuming they're DEL, but DUP/INV/etc also exist..
            new_entry.alts = ['<DEL>']
            new_entry.info["SEGCNT"] = segcnt
            for samp, dat in zip(new_entry.samples, reg.data[1]):
                if dat == 0:
                    new_entry.samples[samp]["GT"] = (0, 0)
                elif dat == 1:
                    new_entry.samples[samp]["GT"] = (0, 1)
                else:
                    new_entry.samples[samp]["GT"] = (1, 1)
                new_entry.samples[samp]["SEG"] = int(dat)

            out.write(new_entry)
