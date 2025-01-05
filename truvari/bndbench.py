"""
BND SV benchmarking tool
Given baseline and comparison sets of variants, calculate the recall/precision/f-measure
"""
import os
import re
import sys
import logging
import argparse

import pysam
import numpy as np

import truvari
from truvari.bench import BenchOutput, pick_single_matches

def parse_args(args):
    """
    Parse list of command line arguments
    """
    parser = argparse.ArgumentParser(prog="bndbench", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-b", "--base", type=str, required=True,
                        help="Baseline truth-set calls")
    parser.add_argument("-c", "--comp", type=str, required=True,
                        help="Comparison set of calls")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output directory")
    parser.add_argument("--includebed", type=str, default=None,
                        help="Bed file of regions to analyze. Only calls within regions are counted")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Verbose logging")
    parser.add_argument("-r", "--refdist", type=truvari.restricted_int, default=100,
                        help="Max reference distance (%(default)s)")
    parser.add_argument("--bSample", type=str, default=None,
                        help="Baseline calls' sample to use (first)")
    parser.add_argument("--cSample", type=str, default=None,
                        help="Comparison calls' sample to use (first)")
    parser.add_argument("--passonly", action="store_true",
                        help="Only consider calls with FILTER == PASS")
    parser.add_argument("--no-ref", default=None, choices=['a', 'b', 'c'],
                        help="Exclude 0/0 or ./. GT calls from all (a), base (b), or comp (c) vcfs (%(default)s)")
    args = parser.parse_args(args)
    args.base = os.path.abspath(args.base)
    args.comp = os.path.abspath(args.comp)
    args.includebed = os.path.abspath(
        args.includebed) if args.includebed else args.includebed
    return args

def check_params(args):
    """
    Checks parameters, True if an error
    """
    return False

# def annotate_entry(entry, match, header):
# def edit_header(my_vcf):

def bnd_direction_strand(bnd: str) -> tuple:
    """
    Parses a BND ALT string to determine its direction and strand.
     ALT  Meaning
    t[p[ piece extending to the right of p is joined after t
    t]p] reverse comp piece extending left of p is joined after t
    ]p]t piece extending to the left of p is joined before t
    [p[t reverse comp piece extending right of p is joined before t
    
    Note that direction of 'left' means the piece is anchored on the left of the breakpoint

    Args:
        bnd (str): The BND ALT string.

    Returns:
        tuple: A tuple containing the direction ("left" or "right") and strand ("direct" or "complement").
    """
    if bnd.startswith('[') or bnd.endswith('['):
        direction = "left"
    elif bnd.startswith(']') or bnd.endswith(']'):
        direction = "right"
    else:
        raise ValueError(f"Invalid BND ALT format: {bnd}")

    # Determine strand based on the position of the base letter
    if bnd[0] not in '[]':  # Base letter is at the start (before brackets)
        strand = "direct"
    elif bnd[-1] not in '[]':  # Base letter is at the end (after brackets)
        strand = "complement"
    else:
        raise ValueError(f"Invalid BND ALT format: {bnd}")

    return direction, strand

def bnd_position(bnd):
    """
    Extracts the chromosome and position from a BND ALT string.
    
    Args:
        bnd (str): The BND ALT string.
        
    Returns:
        tuple: A tuple containing the chromosome (str) and position (int).
    """

    # Regular expression to match the BND format and extract chrom:pos
    match = re.search(r'[\[\]]([^\[\]:]+):(\d+)[\[\]]', bnd)
    if not match:
        raise ValueError(f"Invalid BND ALT format: {bnd}")
    
    chrom = match.group(1)  # Extract the chromosome
    pos = int(match.group(2))  # Extract the position as an integer
    
    return chrom, pos

def build_matches(base_variants, comp_variants, refdist, chunk_id):
    """
    Builds a match matrix for calls
    """
    # All FPs
    if len(base_variants) == 0:
        fps = []
        for cid, c in enumerate(comp_variants):
            ret = truvari.MatchResult()
            ret.comp = c
            ret.matid = ["", f"{chunk_id}.{cid}"]
            fps.append(ret)
            logging.debug("All FP -> %s", ret)
        return fps

    # All FNs
    if len(comp_variants) == 0:
        fns = []
        for bid, b in enumerate(base_variants):
            ret = truvari.MatchResult()
            ret.base = b
            ret.matid = [f"{chunk_id}.{bid}", ""]
            logging.debug("All FN -> %s", ret)
            fns.append(ret)
        return fns

    match_matrix = []
    for bid, base in enumerate(base_variants):
        base_matches = []
        for cid, comp in enumerate(comp_variants):
            mat = truvari.MatchResult()
            mat.base = base
            mat.comp = comp

            mat.matid = [f"{chunk_id}.{bid}", f"{chunk_id}.{cid}"] 
            mat.state = True
            
            bstart = base.pos
            cstart = comp.pos
            mat.st_dist = bstart - cstart
            mat.state &= mat.score < refdist
            
            b_bnd = bnd_direction_strand(base.alts[0])
            c_bnd = bnd_direction_strand(comp.alts[0])
            mat.state &= b_bnd == c_bnd
            
            b_pos2 = bnd_position(base.alts[0])
            c_pos2 = bnd_position(comp.alts[0])
            dist2 = abs(b_pos2[1] - c_pos2[1])
            mat.state &= b_pos2[0] == c_pos2[0]
            mat.state &= dist2 < refdist

            mat.score = abs(mat.st_dist + dist2) / 2
            logging.debug("Made mat -> %s", mat)
            base_matches.append(mat)

        match_matrix.append(base_matches)

    return np.array(match_matrix)

def bnd_filter(variants):
    """
    Only yield variants which are BNDs
    """
    for i in variants:
        if i.alleles_variant_types[1] == 'BND':
            yield i

def bndbench_main(cmdargs):
    """
    Main entry point
    """
    args = parse_args(cmdargs)
    if check_params(args):
        sys.stderr.write("Couldn't run Truvari. Please fix parameters\n")
        sys.exit(100)

    matcher = truvari.Matcher()
    matcher.params.pctseq = 0
    matcher.params.sizemin = 0
    matcher.params.sizefilt = 0
    matcher.params.refdist = args.refdist
    matcher.params.passonly = args.passonly
    matcher.params.bSample = args.bSample
    matcher.params.cSample = args.cSample
    matcher.params.sizemax = sys.maxsize

    #output = BenchOutput(self, self.matcher)
    # beware annotate_entry, and edit_header

    base = pysam.VariantFile(args.base)
    comp = pysam.VariantFile(args.comp)

    region_tree = truvari.build_region_tree(base, comp, args.includebed)
    truvari.merge_region_tree_overlaps(region_tree)

    base_i = bnd_filter(truvari.region_filter(base, region_tree))
    comp_i = bnd_filter(truvari.region_filter(comp, region_tree))
    
    chunks = truvari.chunker(matcher, ('base', base_i), ('comp', comp_i))
    
    tp = 0
    fp = 0
    fn = 0
    for chunk, m_id in chunks:
        matches = build_matches(chunk['base'], chunk['comp'], args.refdist, m_id)
        if isinstance(matches, list):
            for i in matches:
                if i.base:
                    fn += 1
                else:
                    fp += 1
                print('no match', i)
            continue
        picked = pick_single_matches(matches)
        for i in picked:
            if i.state:
                tp += 1
            else:
                if i.base:
                    fn += 1
                else:
                    fp += 1
            print(i)
            #output.write_match(match)
    print(f"TP: {tp}; FP: {fp}, FN: {fn}")
    #output.close_outputs()
    #logging.info("Stats: %s", json.dumps(output.stats_box, indent=4))
    logging.info("Finished bndbench")

