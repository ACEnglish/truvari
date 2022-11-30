"""
Count truvari bench variants per-region
Output bed file with new columns of the variant counts (tpbase, tp, fn, fp)
"""
import os
import sys
import argparse

import pysam
import pandas as pd

import truvari

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="stratify", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("benchdir", metavar="DIR",
                        help="Truvari bench directory")
    parser.add_argument("bedfile", metavar="BED",
                        help="Regions to process")
    parser.add_argument("-o", "--output", metavar="OUT", default="stratifications.bed",
                        help="Output bed-like file")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug, show_version=True)
    return args

def stratify(benchdir, bedfile):
    """
    Count the number of variants per bed region in Truvari bench directory by state
    Returns a pandas dataframe of the bedfile with extra columns of the counts
    """
    regions = pd.read_csv(bedfile, sep='\t', header=None)
    vcfs = {'tpbase': pysam.VariantFile(os.path.join(benchdir, 'tp-base.vcf.gz')),
            'tp': pysam.VariantFile(os.path.join(benchdir, 'tp-call.vcf.gz')),
            'fn': pysam.VariantFile(os.path.join(benchdir, 'fn.vcf.gz')),
            'fp': pysam.VariantFile(os.path.join(benchdir, 'fp.vcf.gz'))
            }
    new_columns = []
    for _, row in regions.iterrows():
        counts = {"tpbase":0, "tp":0, "fn":0, "fp":0}
        for name, fh in vcfs.items():
            m_count = 0
            for entry in fh.fetch(row[0], row[1], row[2]):
                counts[name] += 1
        new_columns.append(counts)
    new_columns = pd.DataFrame(new_columns, index=regions.index)[["tpbase", "tp", "fn", "fp"]]
    regions = regions.join(new_columns)
    return regions

def stratify_main(cmdargs):
    """
    stratify
    """
    args = parse_args(cmdargs)

    regions = stratify(args.benchdir, args.bedfile)
    regions.to_csv(args.output, header=False, index=False, sep='\t')

if __name__ == '__main__':
    stratify_main(sys.argv[1:])
