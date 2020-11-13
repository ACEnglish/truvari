"""
Structural variant collapser

Will collapse all variants that match over thresholds
First variant will be kept in <output.vcf>.
Others will be placed in <output>.collapsed.vcf
Genotypes will be consolidated.
Variants below sizemin will be output directly

I should actually edit the comp entry so that I can give it info on why it matched
"""
# pylint: disable=too-many-statements, no-member
import os
import sys
import json
import logging
import argparse
from collections import defaultdict, namedtuple

import pysam
import pyfaidx

import truvari

MATCHRESULT = namedtuple("matchresult", ("score seq_similarity size_similarity "
                                         "ovl_pct size_diff start_distance "
                                         "end_distance match_entry"))

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
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output vcf")
    parser.add_argument("-f", "--reference", type=str, default=None,
                        help="Indexed fasta used to call variants")
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

    return args

def edit_header(my_vcf):
    """
    Add INFO for new fields to vcf
    #Probably want to put in the PG whatever, too
    """
    # Update header
    # Edit Header
    header = my_vcf.header.copy()
    header.add_line(('##INFO=<ID=NumMerged,Number=1,Type=Integer,'
                     'Description="Number of calls collapsed into this call by truvari">'))
    header.add_line(('##INFO=<ID=CollapseId,Number=1,Type=Integer,'
                     'Description="Truvari uid to help tie output.vcf and output.collapsed.vcf entries together">'))
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
    outputs["seek_vcf"] = pysam.VariantFile(args.input, 'r')
    outputs["header"] = edit_header(outputs["input_vcf"])
    outputs["output_vcf"] = pysam.VariantFile(args.output, 'w', header=outputs["header"])
    merfile = args.output
    if merfile.endswith(".vcf"):
        merfile = merfile[:-4] + ".collapsed.vcf"
    outputs["collap_vcf"] = pysam.VariantFile(merfile, 'w', header=outputs["header"])
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
    new_entry = truvari.copy_entry(entry, outputs["header"])
    new_entry.info["CollapseId"] = match_id
    return new_entry

def edit_output_entry(entry, neighs, match_id, outputs):
    """
    Edit the representative entry that's going to placed in the output vcf
    """
    new_entry = truvari.copy_entry(entry, outputs["header"])
    new_entry.info["CollapseId"] = match_id
    new_entry.info["NumMerged"] = len(neighs)

    # Update with the first genotyped sample's information
    for samp_name in new_entry.samples:
        fmt = new_entry.samples[samp_name]
        m_gt = truvari.stats.get_gt(fmt["GT"])
        if m_gt == truvari.GT.NON:
            idx = 0
            assigned = False
            while not assigned and idx < len(neighs):
                s_fmt = neighs[idx].samples[samp_name]
                s_gt = truvari.stats.get_gt(s_fmt["GT"])
                if s_gt != truvari.GT.NON:
                    assigned = True
                    for key in fmt:
                        fmt[key] = s_fmt[key]
                idx += 1
    return new_entry


def collapse_main(cmdargs):
    """
    Main entry point for running Truvari collapse
    """
    args = parse_args(cmdargs)

    if check_params(args):
        sys.stderr.write("Couldn't run Truvari. Please fix parameters\n")
        sys.exit(100)

    # We can now 'safely' perform everything
    outputs = setup_outputs(args)
    reference = pysam.FastaFile(args.reference) if args.reference else None

    # Calls that have been matched up
    matched_calls = defaultdict(bool)

    # for variant in base - do filtering on it and then try to match it to comp
    logging.info("Merging VCF")
    output_cnt = 0
    kept_cnt = 0
    collap_cnt = 0
    for eid, base_entry in enumerate(outputs["input_vcf"]):
        base_key = truvari.entry_to_key('b', base_entry)
        if base_key in matched_calls:
            continue
        
        # First time seeing it, last time using it
        matched_calls[base_key] = True
        sizeA = truvari.entry_size(base_entry)
        if sizeA < args.sizemin or sizeA > args.sizemax or (args.passonly and truvari.filter_value(base_entry)):
            outputs["output_vcf"].write(base_entry)
            output_cnt += 1
            continue
        
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

            mat = match_calls(base_entry, comp_entry, astart, aend, sizeA, sizeB, reference, args, outputs)
            if isinstance(mat, bool):
                continue
            matched_calls[comp_key] = True
            thresh_neighbors.append(edit_collap_entry(comp_entry, mat, eid, outputs))
        
        if thresh_neighbors:
            new_entry = edit_output_entry(base_entry, thresh_neighbors, eid, outputs)
            outputs["output_vcf"].write(new_entry)
            output_cnt += 1
            kept_cnt += 1
            for i in thresh_neighbors:
                outputs['collap_vcf'].write(i)
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
