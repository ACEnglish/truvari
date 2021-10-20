import sys
import logging
import argparse
import itertools
from functools import cmp_to_key

import pysam
import pandas as pd

import truvari
import truvari.bench as trubench

# reuse chunker and mapper


def collapse_chunk(chunk):
    """
    """
    matcher, chunk_dict, chunk_id = chunk
    logging.debug(f"Comparing chunk {chunk_dict}")
    match_array = []
    for b1, b2 in itertools.combinations(chunk_dict['base'], 2):
        bid1 = chunk_dict['base'].index(b1)
        bid2 = chunk_dict['base'].index(b2)
        mat = matcher.build_match(b1, b2, f'{chunk_id}.{bid1}.{bid2}')
        if matcher.hap:
            if not hap_resolve(b1, b2):
                mat.state = False
        match_array.append(mat)
    match_array.sort(reverse=True, key=matcher.sorter)
    return match_array


def chain_calls():
    """
        given an entry,
        for every hit in the matrix that has this entry as the base
        if its state is True
        Set it to be output to collapse, and do the same chain_calls on this comp.
        Just sets all the bases to None?
    """
    pass


def hap_resolve(entryA, entryB):
    """
    Returns true if the calls' genotypes suggest it can be collapsed
    i.e. if either call is HOM, they can't be collapsed.
    If calls are on the same haplotype (1/0 & 1/0), they cannot be collapsed
    """
    gtA = entryA.samples[0]["GT"]
    gtB = entryB.samples[0]["GT"]
    if gtA == (1, 1) or gtB == (1, 1):
        return False
    if gtA == gtB:
        return False
    return True


def sort_first(b1, b2):
    """
    Sor
    Remove all but the single best neighbor from the list of neighs in-place
    """
    return b1 < b2
    max_score = 0
    max_score_pos = 0
    for pos, n in enumerate(neighs):
        if n[1].score > max_score:
            max_score = n[1].score
            max_score_pos = pos
    return [neighs[max_score_pos]]


def sort_maxqual(b1, b2):
    """
    Swap the entry with the one containing the call with the maximum quality score
    """
    return b1 < b2
    max_qual = entry.qual
    max_qual_idx = -1
    for pos, val in enumerate(neighs):  # pylint: disable=unused-variable
        cmp_qual = neighs[pos].entry.qual
        if cmp_qual > max_qual:
            max_qual_idx = pos

    if max_qual_idx == -1:
        base_entry = entry
    else:
        swap = neighs[max_qual_idx]
        base_entry = swap.entry
        key = truvari.entry_to_key(entry, bounds=True)
        neighs[max_qual_idx] = COLLAPENTRY(
            entry, swap.match, swap.match_id, key)
    return base_entry, neighs


def sort_common(b1, b2):
    """
    Swap the entry with the one containing the highest MAC
    """
    return b1 < b2
    most_mac = truvari.allele_freq_annos(entry)["MAC"]
    most_mac_idx = -1
    for pos, val in enumerate(neighs):  # pylint: disable=unused-variable
        cmp_mac = truvari.allele_freq_annos(neighs[pos].entry)["MAC"]
        if cmp_mac > most_mac:
            most_mac_idx = pos

    if most_mac_idx == -1:
        base_entry = entry
    else:
        swap = neighs[most_mac_idx]
        base_entry = swap.entry
        key = truvari.entry_to_key(entry, bounds=True)
        neighs[most_mac_idx] = COLLAPENTRY(
            entry, swap.match, swap.match_id, key)
    return base_entry, neighs


SORTS = {'first': cmp_to_key(sort_first),
         'maxqual': cmp_to_key(sort_maxqual),
         'common': cmp_to_key(sort_common)}


#######
# VCF #
#######
def edit_header(my_vcf):
    """
    Add INFO for new fields to vcf
    """
    header = my_vcf.header.copy()
    header.add_line(('##INFO=<ID=NumCollapsed,Number=1,Type=Integer,'
                     'Description="Number of calls collapsed into this call by truvari">'))
    header.add_line(('##INFO=<ID=CollapseId,Number=1,Type=Integer,'
                     'Description="Truvari uid to help tie output.vcf and output.collapsed.vcf entries together">'))
    return header

# Collapsed output


def annotate_entry(entry):
    """
    Edit an entry that's going to be collapsed into another entry
    """
    # MatchId
    # NumCollapsed
    new_entry = truvari.copy_entry(entry, outputs["c_header"])
    new_entry.info["CollapseId"] = match_id
    new_entry.info["NumCollapsed"] = match.collapse_count
    return new_entry

##################
# Args & Outputs #
##################


def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="collapse", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Comparison set of calls")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output vcf (stdout)")
    parser.add_argument("-c", "--collapsed-output", type=str, default="collapsed.vcf",
                        help="Where collapsed variants are written (collapsed.vcf)")
    parser.add_argument("-f", "--reference", type=str, default=None,
                        help="Indexed fasta used to call variants")
    parser.add_argument("-k", "--keep", choices=["first", "maxqual", "common"], default="first",
                        help="When collapsing calls, which one to keep (%(default)s)")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Verbose logging")

    # trubench.add_comparison_args(parser)
    thresg = parser.add_argument_group("Comparison Threshold Arguments")
    thresg.add_argument("-r", "--refdist", type=int, default=500,
                        help="Max reference location distance (%(default)s)")
    thresg.add_argument("-p", "--pctsim", type=truvari.restricted_float, default=0.95,
                        help="Min percent allele sequence similarity. Set to 0 to ignore. (%(default)s)")
    thresg.add_argument("-B", "--buffer", type=truvari.restricted_float, default=0.10,
                        help="Percent of the reference span to buffer the haplotype sequence created")
    thresg.add_argument("-P", "--pctsize", type=truvari.restricted_float, default=0.95,
                        help="Min pct allele size similarity (minvarsize/maxvarsize) (%(default)s)")
    thresg.add_argument("-O", "--pctovl", type=truvari.restricted_float, default=0.0,
                        help="Min pct reciprocal overlap (%(default)s) for DEL events")
    thresg.add_argument("-t", "--typeignore", action="store_true", default=False,
                        help="Variant types don't need to match to compare (%(default)s)")
    thresg.add_argument("--use-lev", action="store_true",
                        help="Use the Levenshtein distance ratio instead of edlib editDistance ratio (%(default)s)")

    parser.add_argument("--hap", action="store_true", default=False,
                        help="Collapsing a single individual's haplotype resolved calls (%(default)s)")
    parser.add_argument("--chain", action="store_true", default=False,
                        help="Chain comparisons to extend possible collapsing (%(default)s)")
    parser.add_argument("--null-consolidate", type=str, default=None,
                        help=("Comma separated list of FORMAT fields to consolidate into the kept "
                              "entry by taking the first non-null from all neighbors (%(default)s)"))
    filteg = parser.add_argument_group("Filtering Arguments")
    filteg.add_argument("-s", "--sizemin", type=int, default=50,
                        help="Minimum variant size to consider for comparison (%(default)s)")
    filteg.add_argument("-S", "--sizemax", type=int, default=50000,
                        help="Maximum variant size to consider for comparison (%(default)s)")
    filteg.add_argument("--passonly", action="store_true", default=False,
                        help="Only consider calls with FILTER == PASS")

    args = parser.parse_args(args)
    if args.pctsim != 0 and not args.reference:
        parser.error("--reference is required when --pctsim is set")

    if args.null_consolidate is not None:
        args.null_consolidate = args.null_consolidate.split(',')

    return args


def check_params(args):
    """
    Checks parameters as much as possible.
    All errors are written to stderr without logging since failures mean no output
    """
    check_fail = False
    if not os.path.exists(args.input):
        check_fail = True
        logging.error("File %s does not exist", args.input)
    if not args.input.endswith(".gz"):
        check_fail = True
        logging.error(
            "Input vcf %s does not end with .gz. Must be bgzip'd", args.input)
    if not os.path.exists(args.input + '.tbi'):
        check_fail = True
        logging.error(
            "Input vcf index %s.tbi does not exist. Must be indexed", args.input)
    if args.hap and args.chain:
        check_fail = True
        logging.error("Cannot specify both --hap and --chain")
    if args.hap and args.keep != "first":
        check_fail = True
        logging.error("Using --hap must use --keep first")
    return check_fail


def build_collapse_matcher(args):
    """
    The matcher collapse needs isn't 100% identical to bench
    So set it up for us here
    """
    args.chunksize = args.refdist
    args.gtcomp = False
    args.bSample = None
    args.cSample = None
    args.sizefilt = args.sizemin
    args.no_ref = False
    args.multimatch = False
    matcher = truvari.Matcher(args=args)
    matcher.params.includebed = None
    matcher.keep = args.keep
    matcher.hap = args.hap
    matcher.chain = args.chain
    matcher.sorter = SORTS[args.keep]

    return matcher


def setup_outputs(args):
    """
    Makes all of the output files
    return a ... to get to each of the
    """
    truvari.setup_logging(args.debug)
    logging.info("Params:\n%s", json.dumps(vars(args), indent=4))

    outputs = {}

    outputs["input_vcf"] = pysam.VariantFile(args.input, 'r')
    # These need to be separate, there's only one seek position in a file handler
    outputs["seek_vcf"] = pysam.VariantFile(args.input, 'r')
    outputs["o_header"] = output_edit_header(outputs["input_vcf"])
    outputs["c_header"] = collapse_edit_header(outputs["input_vcf"])
    num_samps = len(outputs["o_header"].samples)
    if args.hap and num_samps != 1:
        logging.error(
            "--hap mode requires exactly one sample. Found %d", num_samps)
        sys.exit(100)
    outputs["output_vcf"] = pysam.VariantFile(
        args.output, 'w', header=outputs["o_header"])
    outputs["collap_vcf"] = pysam.VariantFile(
        args.collapsed_output, 'w', header=outputs["c_header"])
    return outputs


def close_outputs(outputs):
    """
    Close all the files
    """
    outputs["input_vcf"].close()
    outputs["seek_vcf"].close()
    outputs["output_vcf"].close()
    outputs["collap_vcf"].close()


def collapse_main(args):
    args = parse_args(args)
    # Args that Matcher wants but collapser doesn't use
    # check params
    matcher = build_collapse_matcher(args)
    # setup outputs
    base = pysam.VariantFile(args.input)

    chunks = trubench.chunker(matcher, ('base', base))
    for call in itertools.chain.from_iterable(map(collapse_chunk, chunks)):
        print(call)
        # print(len(chunk['base']))
    # parse args
    # setup outputs
    # for call in chunker collapse_calls
    #   output call
    # Done


if __name__ == '__main__':
    collapse_main(sys.argv[1:])
