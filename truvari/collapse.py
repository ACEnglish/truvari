"""
Structural variant collapser

Will collapse all variants within sizemin/max that match over thresholds
All variants outside size boundaries will be placed into the output

When collapsing, the first variant from a matching set of variants will 
be written to the output while the others will be placed in collapsed output.

Samples with no genotype information in the first variant will be filled by the first
collapsed variant containing genotype information.

ToDo:
When using --hap, we assume phased variants from a single individual. Only the
single best non-exact matching call from the other haplotype will be collapsed,
and the consolidated genotype will become 1/1

For example, if we collapse anything at the same position:

    chr1 1 .. GT 0|1
    chr1 1 .. GT 1|0
    chr1 2 .. GT 1|0

will become:
    chr1 1 .. GT 1/1
    chr1 2 .. GT 1|0

When using --chain mode, instead of collapsing all variants matching the first variant
together, we'll collapse all variants in a matching set together. 
For example, if we have

    chr1 5 ..
    chr1 6 ..
    chr1 7 ..

When we collapse anything within 1bp of each other, without --chain, we output:

    chr1 5 ..
    chr1 7 ..

With --chain, we would collapse `chr1 7` as well, producing

    chr1 5 ..

# Just turn this on by default... damn..
When using --detail, we'll record detailed matching information into the 
collapsed-output VCF entrie's infor fields (e.g. PctSim)


--null-consolidate is for special fields that should be consolidated from the FORMAT, for example,
imagine there is a field (FL) that may be null in some entries for a sample, but we want to preserve
the non-nulls in our kept entry:

    chr1 5 .. GT:GQ:FL 0/1:32:.
    chr1 7 .. GT:GQ:FL 0/0:.:FLAG

Without `--null-consolidate FL`, that FLAG will be lost because the call is present in the @5 call (0/1),
so no conslidation of the collapsed entry's SAMPLE information is pulled. But with `--null-consolidate FL`,
even though there is no consolidation triggered by the GT, we'll still preserve the FL:
    
    chr1 5 .. GT:GQ:FL 0/1:32:FLAG

This only works with Number=1/Number=0 FORMAT fields currently
"""
# pylint: disable=too-many-statements, no-member
import os
import sys
import json
import logging
import argparse
from functools import cmp_to_key
from collections import defaultdict, namedtuple

import pysam
import pyfaidx

import truvari
from truvari.bench import edit_header as bench_edit_header

MATCHRESULT = namedtuple("matchresult", ("score seq_similarity size_similarity "
                                         "ovl_pct size_diff start_distance "
                                         "end_distance match_entry"))
COLLAPENTRY = namedtuple("collapentry", ("entry match match_id key"))
   
def parse_args(args):
    """
    Pull the command line parameters
    """
    def restricted_float(x):
        x = float(x)
        if x < 0.0 or x > 1.0:
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
        return x

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
    parser.add_argument("-k", "--keep", choices=["first", "maxqual"], default="first",
                        help="When collapsing calls, which one to keep (%(default)s)")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Verbose logging")

    thresg = parser.add_argument_group("Comparison Threshold Arguments")
    thresg.add_argument("-r", "--refdist", type=int, default=500,
                        help="Max reference location distance (%(default)s)")
    thresg.add_argument("-p", "--pctsim", type=restricted_float, default=0.95,
                        help="Min percent allele sequence similarity. Set to 0 to ignore. (%(default)s)")
    thresg.add_argument("-B", "--buffer", type=restricted_float, default=0.10,
                        help="Percent of the reference span to buffer the haplotype sequence created")
    thresg.add_argument("-P", "--pctsize", type=restricted_float, default=0.95,
                        help="Min pct allele size similarity (minvarsize/maxvarsize) (%(default)s)")
    thresg.add_argument("-O", "--pctovl", type=restricted_float, default=0.0,
                        help="Minimum pct reciprocal overlap (%(default)s) for DEL events")
    thresg.add_argument("-t", "--typeignore", action="store_true", default=False,
                        help="Variant types don't need to match to compare (%(default)s)")
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

def output_edit_header(my_vcf, header=None):
    """
    Add INFO for new fields to vcf
    #Probably want to put in the PG whatever, too
    """
    # Edit Header
    if not header:
        header = my_vcf.header.copy()
    header.add_line(('##INFO=<ID=NumCollapsed,Number=1,Type=Integer,'
                     'Description="Number of calls collapsed into this call by truvari">'))
    header.add_line(('##INFO=<ID=CollapseId,Number=1,Type=Integer,'
                     'Description="Truvari uid to help tie output.vcf and output.collapsed.vcf entries together">'))
    return header

def collapse_edit_header(my_vcf):
    """
    Add output_edit_header and info lines from truvari bench
    """
    header = bench_edit_header(my_vcf)
    header = output_edit_header(my_vcf, header)
    return header

def annotate_tp(entry, match_result):
    """
    Add the matching annotations to a vcf entry
    match_score, match_pctsim, match_pctsize, match_ovlpct, match_szdiff, \
                    match_stdist, match_endist, match_entry
    """
    entry.info["PctSeqSimilarity"] = match_result.seq_similarity
    entry.info["PctSizeSimilarity"] = match_result.size_similarity
    entry.info["PctRecOverlap"] = match_result.ovl_pct
    entry.info["SizeDiff"] = match_result.size_diff
    entry.info["StartDistance"] = match_result.start_distance
    entry.info["EndDistance"] = match_result.end_distance

def check_params(args):
    """
    Checks parameters as much as possible.
    All errors are written to stderr without logging since failures mean no output
    """
    check_fail = False
    if os.path.exists(args.output):
        logging.error("Output '%s' already exists", args.output)
        check_fail = True
    if not os.path.exists(args.input):
        check_fail = True
        logging.error("File %s does not exist", args.input)
    if not args.input.endswith(".gz"):
        check_fail = True
        logging.error("Input vcf %s does not end with .gz. Must be bgzip'd", args.input)
    if not os.path.exists(args.input + '.tbi'):
        check_fail = True
        logging.error("Input vcf index %s.tbi does not exist. Must be indexed", args.input)
    if args.hap and args.chain:
        check_fail = True
        logging.error("Cannot specify both --hap and --chain")
    if args.hap and args.keep != "first":
        check_fail = True
        logging.error("Using --hap must use --keep first")
    return check_fail

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
        logging.error("--hap mode requires exactly one sample. Found %d", num_samps)
        exit(100)
    outputs["output_vcf"] = pysam.VariantFile(args.output, 'w', header=outputs["o_header"])
    outputs["collap_vcf"] = pysam.VariantFile(args.collapsed_output, 'w', header=outputs["c_header"])
    return outputs

def match_calls(base_entry, comp_entry, astart, aend, sizeA, sizeB, reference, args, outputs):
    """
    Compare the base and comp entries.
    We provied astart...sizeA because we've presumably calculated it before
    Note - This is the crucial component of matching.. so needs to be better
    pulled apart for reusability and put into comparisons
    """
    bstart, bend = truvari.entry_boundaries(comp_entry)
    if not truvari.overlaps(astart - args.refdist, aend + args.refdist, bstart, bend):
        return False

    # Someone in the Base call's neighborhood, we'll see if it passes comparisons
    if not args.typeignore and not truvari.same_variant_type(base_entry, comp_entry):
        logging.debug("%s and %s are not the same SVTYPE", str(base_entry), str(comp_entry))
        return True

    size_similarity, size_diff = truvari.sizesim(sizeA, sizeB)
    if size_similarity < args.pctsize:
        logging.debug("%s and %s size similarity is too low (%f)", str(base_entry),
                      str(comp_entry), size_similarity)
        return True

    ovl_pct = truvari.reciprocal_overlap(astart, aend, bstart, bend)
    if truvari.entry_variant_type(base_entry) == "DEL" and ovl_pct < args.pctovl:
        logging.debug("%s and %s overlap percent is too low (%f)", str(base_entry), str(comp_entry), ovl_pct)
        return True

    if args.pctsim > 0:
        seq_similarity = truvari.entry_pctsim_lev(base_entry, comp_entry, reference, buf_len=args.buffer)
        if seq_similarity < args.pctsim:
            logging.debug("%s and %s sequence similarity is too low (%f)", str(
                base_entry), str(comp_entry), seq_similarity)
            return True
    else:
        seq_similarity = 0

    start_distance = astart - bstart
    end_distance = aend - bend

    score = truvari.weighted_score(seq_similarity, size_similarity, ovl_pct)

    return MATCHRESULT(score, seq_similarity, size_similarity, ovl_pct, size_diff,
                       start_distance, end_distance, comp_entry)

def close_outputs(outputs):
    """
    Close all the files
    """
    outputs["input_vcf"].close()
    outputs["seek_vcf"].close()
    outputs["output_vcf"].close()
    outputs["collap_vcf"].close()

def edit_collap_entry(entry, match, match_id, outputs):
    """
    Edit an entry that's going to be collapsed into another entry
    """
    new_entry = truvari.copy_entry(entry, outputs["c_header"])
    new_entry.info["CollapseId"] = match_id
    new_entry.info["TruScore"] = match.score
    truvari.bench.annotate_tp(new_entry, match)
    return new_entry

def same_hap(entryA, entryB):
    """
    Returns if two entry's first sample is on the same haplotype
    by using GT phasing information
    """
    gtA = entryA.samples[0]["GT"]
    gtB = entryB.samples[0]["GT"]
    return gtA == gtB

def select_best(neighs):
    """
    Remove all but the single best neighbor from the list of neighs in-place
    """
    max_score = 0
    max_score_pos = 0
    for pos, n in enumerate(neighs):
        if n[1].score > max_score:
            max_score = n[1].score
            max_score_pos = pos
    i = len(neighs) - 1
    while i > max_score_pos:
        neighs.pop(i)
        i -= 1
    for i in range(max_score_pos):
        neighs.pop(0)

def edit_output_entry(entry, neighs, match_id, hap, outputs, nullconso=None):
    """
    Edit the representative entry that's going to placed in the output vcf
    If hap, consolidate with haplotype mode and edit the neighs in-place 
    down to only the single best match.
    Force is a set of FORMAT fields that will be checked for nulls and will overwrite between entries
    """
    if nullconso is None:
        nullconso = []

    new_entry = truvari.copy_entry(entry, outputs["o_header"])
    new_entry.info["CollapseId"] = match_id
    new_entry.info["NumCollapsed"] = len(neighs)

    if hap:
        select_best(neighs)
        new_entry.samples[0]["GT"] = (1, 1)
        return new_entry


    # Update with the first genotyped sample's information
    for samp_name in new_entry.samples:
        fmt = new_entry.samples[samp_name]
        m_gt = truvari.stats.get_gt(fmt["GT"])
        # Update the null_consolidates first - these will replace with the first non-null
        for key in nullconso:
            if fmt[key] is None:
                idx = 0
                assigned = False
                while not assigned and idx < len(neighs):
                    s_fmt = neighs[idx][0].samples[samp_name]
                    if key in s_fmt and s_fmt[key] is not None:
                        assigned = True
                        fmt[key] = s_fmt[key]
                    idx += 1
        
        # Now pull the rest based on GT presence
        if m_gt in [truvari.GT.NON, truvari.GT.REF]:
            idx = 0
            assigned = False
            while not assigned and idx < len(neighs):
                s_fmt = neighs[idx][0].samples[samp_name]
                s_gt = truvari.stats.get_gt(s_fmt["GT"])
                if s_gt not in [truvari.GT.NON, truvari.GT.REF]:
                    assigned = True
                    for key in [_ for _ in fmt if _ not in nullconso]:
                        try:
                            fmt[key] = s_fmt[key]
                        except Exception as e:
                            pass
                idx += 1
    return new_entry

def find_neighbors(base_entry, match_id, reference, matched_calls, args, outputs):
    """
    Find all matching neighbors of a call
    returns a list of (edited_neighbors, match_id, 
    """
    base_entry_size = truvari.entry_size(base_entry)
    thresh_neighbors = []
    astart, aend = truvari.entry_boundaries(base_entry)
    fetch_start = max(0, astart - args.refdist - 1)
    fetch_end = aend + args.refdist + 1
    for comp_entry in outputs['seek_vcf'].fetch(base_entry.chrom, fetch_start, fetch_end):
        comp_key = truvari.entry_to_key('b', comp_entry)
        # Don't collapse with anything that's already been seen
        # Which will include looking at itself
        if matched_calls[comp_key]:
            continue
        
        sizeB = truvari.entry_size(comp_entry)
        if sizeB < args.sizemin or sizeB > args.sizemax or (args.passonly and truvari.filter_value(comp_entry)):
            continue

        if args.hap and same_hap(base_entry, comp_entry):
            continue

        mat = match_calls(base_entry, comp_entry, astart, aend, base_entry_size, sizeB, reference, args, outputs)
        if isinstance(mat, bool):
            continue

        # When not in haplotype mode
        # mark this entry as being matched already
        if not args.hap:
            matched_calls[comp_key] = True
        
        thresh_neighbors.append(COLLAPENTRY(comp_entry, mat, match_id, comp_key))
    return thresh_neighbors

def select_maxqual(entry, neighs):
    """
    Swap the entry with the one containing the call with the maximum quality score
    """
    max_qual = entry.qual
    max_qual_idx = -1
    for pos, i in enumerate(neighs):
        cmp_qual = neighs[pos].entry.qual
        if cmp_qual > max_qual:
            max_qual_idx = pos
    
    if max_qual_idx == -1:
        base_entry = entry
    else:
        swap = neighs[max_qual_idx]
        base_entry = swap.entry
        key = truvari.entry_to_key('b', entry)
        neighs[max_qual_idx] = COLLAPENTRY(entry, swap.match, swap.match_id, key)
    return base_entry, neighs

def collapse_main(cmdargs):
    """
    Main entry point for running Truvari collapse
    """
    args = parse_args(cmdargs)

    if check_params(args):
        sys.stderr.write("Couldn't run Truvari. Please fix parameters\n")
        sys.exit(100)

    outputs = setup_outputs(args)
    reference = pysam.FastaFile(args.reference) if args.reference else None

    # Calls that have been matched up
    matched_calls = defaultdict(bool)

    # for variant in base - do filtering on it and then try to match it to comp
    logging.info("Collapsing VCF")
    output_cnt = 0
    kept_cnt = 0
    collap_cnt = 0
    for match_id, base_entry in enumerate(outputs["input_vcf"]):
        base_key = truvari.entry_to_key('b', base_entry)
        if matched_calls[base_key]:
            continue
        
        sizeA = truvari.entry_size(base_entry)
        if sizeA < args.sizemin or sizeA > args.sizemax or (args.passonly and truvari.filter_value(base_entry)):
            outputs["output_vcf"].write(base_entry)
            output_cnt += 1
            continue

        # First time seeing it, last time using it
        matched_calls[base_key] = True
        
        thresh_neighbors = []
        new_neighs = find_neighbors(base_entry, match_id, reference, matched_calls, args, outputs)
        thresh_neighbors.extend(new_neighs)
        # For all the neighbors, we want to go also collapse anything matching them
        while args.chain and new_neighs:
            next_new_neighs = []
            for i in new_neighs:
                next_new_neighs.extend(find_neighbors(i[0], match_id, reference, matched_calls, args, outputs))
            thresh_neighbors.extend(next_new_neighs)
            new_neighs = next_new_neighs

        if args.keep == "maxqual":
            base_entry, thresh_neighbors = select_maxqual(base_entry, thresh_neighbors)
            
        if thresh_neighbors:
            new_entry = edit_output_entry(base_entry, thresh_neighbors, match_id, args.hap, outputs, args.null_consolidate)
            outputs["output_vcf"].write(new_entry)
            
            # Also want to record the original beside its' collapse neighbors
            new_entry = truvari.copy_entry(base_entry, outputs["c_header"])
            new_entry.info["CollapseId"] = match_id
            new_entry.info["NumCollapsed"] = len(thresh_neighbors)
            outputs["collap_vcf"].write(new_entry)
            
            output_cnt += 1
            kept_cnt += 1
            for neigh in thresh_neighbors:
                matched_calls[neigh.key] = True
                out_entry = edit_collap_entry(neigh.entry, neigh.match, neigh.match_id, outputs)
                outputs['collap_vcf'].write(out_entry)
            collap_cnt += len(thresh_neighbors)
        else:
            outputs["output_vcf"].write(base_entry)
            output_cnt += 1

        # Finished with this base entry
    logging.info("%d calls written to output", output_cnt)
    logging.info("%d collapsed calls", collap_cnt)
    logging.info("%d representative calls kept", kept_cnt)
    
    # Close to flush vcfs
    close_outputs(outputs)
    logging.info("Finished")

if __name__ == '__main__':
    collapse_main(sys.argv[1:])
