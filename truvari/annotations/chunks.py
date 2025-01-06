"""
Count the number of variants in each chunk
Column 3: total number of variants
Column 4: comma-deliminted number of sub-chunks after accounting for size
Column 5: comma-deliminted number of sub-chunks after accounting for size and distance again
"""
import sys
import argparse

import truvari
from truvari.collapse import tree_size_chunker, tree_dist_chunker

def parse_args(args):
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(prog="chunks", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", type=str,
                        help="Input VCF")
    parser.add_argument("-b", "--bed", type=str, default=None,
                        help="Bed file of variants to chunk")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output name")
    parser.add_argument("-c", "--chunksize", type=int, default=500,
                        help="Distance between variants to split chunks (%(default)s)")
    parser.add_argument("-s", "--sizemin", type=int, default=50,
                        help="Minimum SV length (%(default)s)")
    parser.add_argument("-S", "--sizemax", type=int, default=50000,
                        help="Maximum SV length (%(default)s)")
    args = parser.parse_args(args)
    truvari.setup_logging(show_version=True)
    return args

def get_bounds(cnk):
    """
    Min start and max end of variants
    """
    mstart = sys.maxsize
    mend = 0
    for i in cnk:
        mstart = min(mstart, i.start)
        mend = max(mend, i.end)
    return mstart, mend

def chunks_main(args):
    """
    Main
    """
    args = parse_args(args)
    v = truvari.VariantFile(args.input)
    m = truvari.Matcher()
    m.params.pctseq = 0
    m.params.sizemin = args.sizemin
    m.params.sizefilt = args.sizemin
    m.params.sizemax = args.sizemax
    m.params.chunksize = args.chunksize
    m.params.refdist = args.chunksize
    if args.bed:
        regions = truvari.build_region_tree(v, includebed=args.bed)
        v = truvari.region_filter(v, regions)
    c = truvari.chunker(m, ('base', v))

    with open(args.output, 'w') as fout:
        for chunk, _ in c:
            if not chunk['base']:
                continue
            s, e = get_bounds(chunk['base'])
            chrom = chunk['base'][0].chrom
            num = len(chunk['base'])
            fout.write(f"{chrom}\t{s}\t{e}\t{num}")
            s_cnts = []
            d_cnts = []
            for i, _ in tree_size_chunker(m, [(chunk, 0)]):
                if i['base']:
                    s_cnts.append(len(i['base']))
                for j, _ in tree_dist_chunker(m, [(i, 0)]):
                    d_cnts.append(len(j['base']))
            s_cnts = ",".join(map(str, s_cnts))
            d_cnts = ",".join(map(str, d_cnts))
            fout.write(f"\t{s_cnts}\t{d_cnts}\n")
