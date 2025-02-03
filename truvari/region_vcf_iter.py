"""
Helper class to specify included regions of the genome when iterating events.
"""
import os
import sys
import copy
import logging
from collections import defaultdict, deque

from intervaltree import IntervalTree
import truvari


def build_region_tree(vcfA, vcfB=None, includebed=None):
    """
    Builds a dictionary of chromosome names mapped to `IntervalTree` objects
    containing regions based on the input VCF files and an optional BED file.

    Parameters
    ----------
    vcfA : pysam.VariantFile
        The primary VCF file used to define the regions and contigs.
    vcfB : pysam.VariantFile, optional
        The secondary VCF file used for comparison. If provided, its contigs
        are considered to exclude any contigs not present in `vcfA`.
        Default is None.
    includebed : str, optional
        Path to a BED file specifying regions to include. If provided,
        these regions override those defined by the VCF files. Contigs in the
        BED file not found in `vcfA` are excluded. Default is None.

    Returns
    -------
    dict
        A dictionary where keys are chromosome names (str), and values are
        `IntervalTree` objects representing the regions.

    Notes
    -----
    - If `includebed` is provided:
        - The function reads regions from the BED file and uses them to define the
          regions in `all_regions`.
        - Contigs in the BED file that are not present in `vcfA` are excluded,
          and a warning is logged.
    - Contigs in `vcfB` that are not present in `vcfA` are ignored, and a
      warning is logged.
    - If a contig in `vcfA` lacks a length definition, the program logs an
      error and exits with code 10.

    Raises
    ------
    SystemExit
        If a contig in `vcfA` has no length defined, the program terminates
        with an error message and exit code 10.

    Logging
    -------
    - Logs warnings if contigs in `vcfB` are not found in `vcfA`.
    - Logs warnings if contigs in the BED file are not found in `vcfA`.
    - Logs an error and exits if any contig in `vcfA` lacks a length definition.

    """
    contigA_set = set(vcfA.header.contigs.keys())
    contigB_set = set(vcfB.header.contigs.keys()) if vcfB else contigA_set

    all_regions = defaultdict(IntervalTree)
    excluding = contigB_set - contigA_set
    if excluding:
        logging.warning(
            "%d contigs present in comparison VCF header are not in baseline VCF.", len(excluding))

    if includebed is not None:
        all_regions, counter = read_bed_tree(includebed)
        logging.info("Including %d bed regions", counter)
        excluding = set(all_regions.keys()) - contigA_set
        if excluding:
            logging.warning(
                "Excluding %d contigs in BED not in baseline VCF.", len(excluding))
        for i in excluding:
            del all_regions[i]
        return all_regions

    for contig in contigA_set:
        name = vcfA.header.contigs[contig].name
        length = vcfA.header.contigs[contig].length
        if not length:
            logging.error(
                "Contig %s has no length definition. Fix header.", name)
            sys.exit(10)
        all_regions[name].addi(0, length + 1)
    return all_regions


def merge_region_tree_overlaps(tree):
    """
    Runs IntervalTree.merge_overlaps on all trees. Info Logs the count of chromosomes with overlapping regions 
    and Debug Logs the chromosome names with overlapping regions.
    """
    chr_with_overlaps = []
    for i in tree:
        pre_len = len(tree[i])
        tree[i].merge_overlaps()
        post_len = len(tree[i])
        if pre_len != post_len:
            chr_with_overlaps.append(i)
    if chr_with_overlaps:
        logging.info("Found %d chromosomes with overlapping regions",
                     len(chr_with_overlaps))
        logging.debug("CHRs: %s", chr_with_overlaps)


def extend_region_tree(tree, pad):
    """
    Extends all intervals by a fixed number of bases on each side
    Returns a copy of this tree
    """
    logging.info("Extending the regions by %d bases", pad)
    n_tree = copy.deepcopy(tree)

    for chrom in n_tree:
        n_tree[chrom] = IntervalTree.from_tuples(
            ((max(0, i.begin - pad), i.end + pad)) for i in n_tree[chrom])

    truvari.merge_region_tree_overlaps(n_tree)
    return n_tree


def read_bed_tree(filename, chrom_col=0, start_col=1, end_col=2, one_based=False, comment='#', idxfmt=None):
    """
    Build an dictionary of IntervalTrees for each chromosome from tab-delimited annotation file

    By default, the file is assumed to be a bed-format. If custom chrom/start/end are used, the columns can be
    specified.

    idxfmt is a string that will be formatted with the line number from the file (e.g. "num %d" will
    make intervals with data="num 0", "num 1". By default the data will be interger line number.
    If intervals will be compared between anno_trees, set idxfmt to ""

    :param `filename`: Path to file to parse, can be compressed
    :type `filename`: string
    :param `chrom_col`: Index of column in file with chromosome
    :type `chrom_col`: integer, optional
    :param `start_col`: Index of column in file with start
    :type `start_col`: integer, optional
    :param `end_col`: Index of column in file with end
    :type `end_col`: integer, optional
    :param `one_based`: True if coordinates are one-based
    :type `one_based`: bool, optional
    :param `is_pyintv`: add 1 to end position to correct pyintervaltree.overlap behavior
    :type `is_pyintv`: bool, optional
    :param `comment`: ignore lines if they start with this string
    :type `comment`: string, optional
    :param `idxfmt`: Index of column in file with chromosome
    :type `idxfmt`: string, optional

    :return: dictionary with chromosome keys and :class:`intervaltree.IntervalTree` values, number of lines parsed
    :rtype: tuple (dict, int)
    """
    idx = 0
    correction = 1 if one_based else 0
    ttree = defaultdict(list)
    for line in truvari.opt_gz_open(filename):
        if line.startswith(comment):
            continue
        data = line.strip().split('\t')
        chrom = data[chrom_col]
        start = int(data[start_col]) - correction
        end = int(data[end_col])
        if idxfmt is not None:
            m_idx = idxfmt.format(idx)
        else:
            m_idx = idx
        ttree[chrom].append((start, end + 1, m_idx))
        idx += 1
    tree = {}
    for chrom in ttree:
        tree[chrom] = IntervalTree.from_tuples(ttree[chrom])
    return tree, idx


def region_filter(vcf, tree, inside=True, with_region=False):
    """
    Chooses to stream or fetch entries inside/outside a VCF
    If the VCF is over 25Mb or the number of regions is above 1k, use stream
    """
    sz = os.stat(vcf.filename)
    if not inside or sz.st_size > (25 * 2**20) or sum(len(_) for _ in tree.values()) > 1000:
        return region_filter_stream(vcf, tree, inside, with_region)

    return region_filter_fetch(vcf, tree, with_region)


def region_filter_fetch(vcf, tree, with_region=False):
    """
    Given a VariantRecord iter and defaultdict(IntervalTree),
    yield variants which are inside/outside the tree regions
    The region associated with the entry can be retuned also when using with_region.
    with_region returns (entry, (chrom, Interval))
    Can only check for variants within a region
    """
    ret_type = (lambda x, y, z: (x, (y, z))
                ) if with_region else (lambda x, y, z: x)
    for chrom in sorted(tree.keys()):
        for intv in sorted(tree[chrom]):
            try:
                for entry in vcf.fetch(chrom, intv.begin, intv.end):
                    if entry.within(intv.begin, intv.end - 1):
                        yield ret_type(entry, chrom, intv)
            except ValueError:
                logging.warning("Unable to fetch %s from %s",
                                chrom, vcf.filename)


def region_filter_stream(vcf, tree, inside=True, with_region=False):
    """
    Given a VariantRecord iter and defaultdict(IntervalTree),
    yield variants which are inside/outside the tree regions
    The region associated with the entry can be retuned also when using with_region.
    with_region returns (entry, (chrom, Interval))
    """
    if not inside:
        chroms = sorted(list(vcf.header.contigs))
    else:
        chroms = sorted(tree.keys())
    ret_type = (lambda x, y, z: (x, (y, z))
                ) if with_region else (lambda x, y, z: x)
    for chrom in chroms:
        cur_tree = deque(sorted(tree[chrom]))
        try:
            cur_intv = cur_tree.popleft()
        except (IndexError, KeyError):
            # region-less chromosome
            if not inside:
                try:
                    for cur_entry in vcf.fetch(chrom):
                        yield ret_type(cur_entry, chrom, cur_intv)
                except ValueError:
                    pass  # region on chromosome not in vcf
            continue

        try:
            cur_iter = vcf.fetch(chrom)
        except ValueError:
            continue  # region on chromosome not in vcf
        try:
            cur_entry = next(cur_iter)
        except StopIteration:
            # variant-less chromosome
            continue
        cur_start, cur_end = cur_entry.boundaries()

        while True:
            # if start is after this interval, we need the next interval
            if cur_start > cur_intv.end:
                try:
                    cur_intv = cur_tree.popleft()
                except IndexError:
                    if not inside:
                        # pass this back before flush after the while
                        yield ret_type(cur_entry, chrom, cur_intv)
                    break
            # well before, we need the next entry
            elif cur_end < cur_intv.begin:
                if not inside:
                    yield ret_type(cur_entry, chrom, cur_intv)
                try:
                    cur_entry = next(cur_iter)
                    cur_start, cur_end = cur_entry.boundaries()
                except StopIteration:
                    break
            else:
                end_within = cur_entry.var_type() != truvari.SV.INS
                is_within = truvari.coords_within(
                    cur_start, cur_end, cur_intv.begin, cur_intv.end - 1, end_within)
                if is_within == inside:
                    yield ret_type(cur_entry, chrom, cur_intv)
                try:
                    cur_entry = next(cur_iter)
                    cur_start, cur_end = cur_entry.boundaries()
                except StopIteration:
                    break

        # if we finished the intervals first, need to flush the rest of the outside entries
        if not inside:
            for cur_entry in cur_iter:
                yield ret_type(cur_entry, chrom, cur_intv)
