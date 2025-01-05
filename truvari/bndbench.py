"""
BND SV benchmarking tool
Given baseline and comparison sets of variants, calculate the recall/precision/f-measure
"""
import re
import logging
import numpy as np
import truvari

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

def build_matches(base_variants, comp_variants, matcher, chunk_id):
    """
    Builds a match matrix for calls
    """
    # All FPs
    if len(base_variants) == 0:
        fps = []
        for cid, comp in enumerate(comp_variants):
            ret = truvari.MatchResult()
            ret.comp = comp
            ret.matid = ["", f"{chunk_id}.{cid}"]
            fps.append(ret)
            logging.debug("All FP -> %s", ret)
        return fps

    # All FNs
    if len(comp_variants) == 0:
        fns = []
        for bid, base in enumerate(base_variants):
            ret = truvari.MatchResult()
            ret.base = base
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
            mat.state = base.chrom == comp.chrom

            mat.st_dist = base.pos - comp.pos
            mat.state &= abs(mat.st_dist) < matcher.params.bnddist
            b_bnd = bnd_direction_strand(base.alts[0])
            c_bnd = bnd_direction_strand(comp.alts[0])
            mat.state &= b_bnd == c_bnd

            b_pos2 = bnd_position(base.alts[0])
            c_pos2 = bnd_position(comp.alts[0])
            mat.ed_dist = b_pos2[1] - c_pos2[1]
            mat.state &= b_pos2[0] == c_pos2[0]
            mat.state &= mat.ed_dist < matcher.params.bnddist

            # not dry
            if "GT" in base.samples[matcher.params.bSample]:
                mat.base_gt = base.samples[matcher.params.bSample]["GT"]
                mat.base_gt_count = sum(1 for _ in mat.base_gt if _ == 1)
            if "GT" in comp.samples[matcher.params.cSample]:
                mat.comp_gt = comp.samples[matcher.params.cSample]["GT"]
                mat.comp_gt_count = sum(1 for _ in mat.comp_gt if _ == 1)
            mat.gt_match = abs(mat.base_gt_count - mat.comp_gt_count)

            # Score is percent of allowed distance needed to find this match
            mat.score = (1 - ((abs(mat.st_dist) + abs(mat.ed_dist)) / 2) / matcher.params.bnddist) * 100
            # I think I'm missing GT stuff here
            logging.debug("BND match -> %s", mat)
            base_matches.append(mat)

        match_matrix.append(base_matches)

    return np.array(match_matrix)
