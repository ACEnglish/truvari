"""
Structural variant collapser

Will collapse all variants within sizemin/max that match over thresholds
"""

import os
import sys
import json
import logging
import argparse
import itertools
from functools import cmp_to_key

import pysam

import truvari
import truvari.bench as trubench


def collapse_chunk(chunk):
    """
    For calls in a chunk, separate which ones should be kept from collapsed
    """
    matcher, chunk_dict, chunk_id = chunk
    calls = chunk_dict['base']
    logging.debug(f"Comparing chunk {calls}")
    calls.sort(reverse=True, key=matcher.sorter)

    # keep_key : [keep entry, [collap entries], match_id]
    ret = {}
    # collap_key : keep_key
    chain_lookup = {}

    call_id = 0
    while calls:
        # Take the first variant
        cur_keep_candidate = calls.pop(0)
        keep_key = truvari.entry_to_key(cur_keep_candidate)
        # if chain, any matches we find will be added to a previous keep
        if matcher.chain and keep_key in chain_lookup:
            keep_key = chain_lookup[keep_key]
        else:  # otherwise, this is assumed to be a keep call
            ret[keep_key] = [cur_keep_candidate, [], f'{chunk_id}.{call_id}']

        # Separate calls collapsing with this keep from the rest
        remaining_calls = []
        for cur_collapse_candidate in calls:
            mat = matcher.build_match(cur_keep_candidate,
                                      cur_collapse_candidate,
                                      ret[keep_key][2])
            if matcher.hap and not hap_resolve(cur_keep_candidate, cur_collapse_candidate):
                mat.state = False
            # to collapse
            if mat.state:
                if matcher.chain:
                    # need to record for downstream
                    collap_key = truvari.entry_to_key(cur_collapse_candidate)
                    chain_lookup[collap_key] = keep_key
                    remaining_calls.append(cur_collapse_candidate)
                ret[keep_key][1].append(mat)
            else:
                remaining_calls.append(cur_collapse_candidate)

        # Only put in the single best match
        # Leave the others to be handled later
        if matcher.hap and ret[keep_key][1]:
            candidates = sorted(ret[keep_key][1], reverse=True)
            ret[keep_key][1] = [candidates.pop(0)]
            remaining_calls.extend(candidates)

        calls = remaining_calls

    for key, val in ret.items():
        logging.debug("Collapsing %s", key)
        val[0] = collapse_into_entry(val[0], val[1])

    ret = list(ret.values())
    for i in chunk_dict['__filtered']:
        ret.append([i, None, None])

    return ret


def collapse_into_entry(entry, others, hap_mode=False):
    """
    Consolidate information for genotypes where sample is unset
    okay - this works, but its a mess
    """
    # short circuit
    if not others:
        return entry

    # We'll populate with the most similar, first
    others.sort(reverse=True)
    # I have a special case of --hap. I need to allow hets
    replace_gts = ["UNK", "REF", "NON"]
    if hap_mode:
        replace_gts.append("HET")

    # Each sample of this entry needs to be checked/set
    for sample in entry.samples:
        m_gt = truvari.get_gt(entry.samples[sample]["GT"]).name
        if m_gt not in replace_gts:
            continue  # already set
        n_idx = None
        for pos, o_entry in enumerate(others):
            o_entry = o_entry.comp
            o_gt = truvari.get_gt(o_entry.samples[sample]["GT"]).name
            if o_gt not in replace_gts:
                n_idx = pos
                break  # this is the first other that's set
        # consolidate
        if hap_mode and m_gt == "HET":
            entry.samples[sample]["GT"] = (1, 1)
        elif n_idx is not None:
            o_entry = others[n_idx].comp
            for key in set(entry.samples[sample].keys() + o_entry.samples[sample].keys()):
                entry.samples[sample][key] = o_entry.samples[sample][key]

    return entry


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
    Remove all but the single best neighbor from the list of neighs in-place
    """
    return b1.pos < b2.pos


def sort_maxqual(b1, b2):
    """
    Swap the entry with the one containing the call with the maximum quality score
    """
    return b1.qual < b2.qual


def sort_common(b1, b2):
    """
    Swap the entry with the one containing the highest MAC
    """
    mac1 = truvari.allele_freq_annos(b1)["MAC"]
    mac2 = truvari.allele_freq_annos(b2)["MAC"]
    return mac1 < mac2


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
    header.add_line(('##INFO=<ID=CollapseId,Number=1,Type=String,'
                     'Description="Truvari uid to help tie output.vcf and output.collapsed.vcf entries together">'))
    return header


def annotate_entry(entry, num_collapsed, match_id, header):
    """
    Edit an entry that's going to be collapsed into another entry
    """
    new_entry = truvari.copy_entry(entry, header)
    new_entry.info["NumCollapsed"] = num_collapsed
    new_entry.info["CollapseId"] = match_id
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

    in_vcf = pysam.VariantFile(args.input)
    outputs["o_header"] = edit_header(in_vcf)
    outputs["c_header"] = trubench.edit_header(in_vcf)
    num_samps = len(outputs["o_header"].samples)
    if args.hap and num_samps != 1:
        logging.error(
            "--hap mode requires exactly one sample. Found %d", num_samps)
        sys.exit(100)
    outputs["output_vcf"] = pysam.VariantFile(args.output, 'w',
                                              header=outputs["o_header"])
    outputs["collap_vcf"] = pysam.VariantFile(args.collapsed_output, 'w',
                                              header=outputs["c_header"])
    outputs["stats_box"] = {"collap_cnt": 0, "kept_cnt": 0, "out_cnt": 0}
    return outputs


def output_writer(to_collapse, outputs):
    """
    Annotate and write kept/collapsed calls to appropriate files
    """
    entry, collapsed, match_id = to_collapse
    # Save copying time
    if not collapsed:
        outputs["output_vcf"].write(entry)
        outputs["stats_box"]["out_cnt"] += 1
        return

    entry = annotate_entry(entry, len(collapsed),
                           match_id, outputs["o_header"])
    outputs["output_vcf"].write(entry)
    outputs["stats_box"]["out_cnt"] += 1
    outputs["stats_box"]["kept_cnt"] += 1
    for match in collapsed:
        entry = trubench.annotate_entry(match.comp, match, outputs["c_header"])
        outputs["collap_vcf"].write(entry)
        outputs['stats_box']["collap_cnt"] += 1


def close_outputs(outputs):
    """
    Close all the files
    """
    outputs["output_vcf"].close()
    outputs["collap_vcf"].close()


def collapse_main(args):
    """
    Main
    """
    args = parse_args(args)

    if check_params(args):
        sys.stderr.write("Couldn't run Truvari. Please fix parameters\n")
        sys.exit(100)

    matcher = build_collapse_matcher(args)
    outputs = setup_outputs(args)
    base = pysam.VariantFile(args.input)

    chunks = trubench.chunker(matcher, ('base', base))
    for call in itertools.chain.from_iterable(map(collapse_chunk, chunks)):
        output_writer(call, outputs)

    close_outputs(outputs)
    logging.info("Wrote %d Variants", outputs["stats_box"]["out_cnt"])
    logging.info("%d variants collapsed into %d variants", outputs["stats_box"]["collap_cnt"],
                 outputs["stats_box"]["kept_cnt"])
    logging.info("Finished collapse")


if __name__ == '__main__':
    collapse_main(sys.argv[1:])
