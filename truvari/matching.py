"""
Comparison engine
"""
import logging
from collections import Counter, defaultdict
from functools import total_ordering
import pysam


@total_ordering
class MatchResult():  # pylint: disable=too-many-instance-attributes
    """
    Holds results from a matching operation

    Attributes
    ----------
    .. list-table::
       :header-rows: 1

       * - Attribute
         - Description
       * - `base`
         - The base (a.k.a. self) variant
       * - `comp`
         - The comp (a.k.a. other) variant
       * - `state`
         - Boolean of if variants match
       * - `seqsim`
         - Sequence similarity of variants
       * - `sizesim`
         - Size similarity of variants
       * - `ovlpct`
         - Reciprocal overlap ov variants
       * - `sizediff`
         - Base variant var_size minus comp variant var_size
       * - `st_dist`
         - Base variant start position minus comp variant start position
       * - `ed_dist`
         - Base variant end position minus comp variant end position
       * - `gt_match`
         -  Boolean of if genotypes match
       * - `score`
         - TruScore of matches
       * - `matid`
         - Place to put MatchIds, not populated by `truvari.VariantRecord.match`

    """
    __slots__ = ["base", "comp", "base_gt", "base_gt_count", "comp_gt", "comp_gt_count",
                 "state", "seqsim", "sizesim", "ovlpct", "sizediff", "st_dist", "ed_dist",
                 "gt_match", "multi", "score", "matid"]

    def __init__(self):
        self.base = None
        self.comp = None
        self.base_gt = None
        self.base_gt_count = 0
        self.comp_gt = None
        self.comp_gt_count = 0
        self.matid = None
        self.seqsim = None
        self.sizesim = None
        self.ovlpct = None
        self.sizediff = None
        self.st_dist = None
        self.ed_dist = None
        self.gt_match = None
        self.multi = None
        self.state = False
        self.score = None

    def calc_score(self):
        """
        Unite the similarity measures and make a score
        """
        if None not in [self.seqsim, self.sizesim, self.ovlpct]:
            self.score = (self.seqsim + self.sizesim + self.ovlpct) / 3.0 * 100

    def __lt__(self, other):
        def s_abs(value):
            """
            Negative because we want them to be closer
            """
            return -abs(value) if value is not None else -float('inf')
        return (
            (self.state, self.score if self.score is not None else 0, s_abs(self.st_dist), s_abs(self.ed_dist)) <
            (other.state, other.score if other.score is not None else 0, s_abs(other.st_dist), s_abs(other.ed_dist))
        )

    def __eq__(self, other):
        return self.state == other.state and self.score == other.score

    def __str__(self):
        return f'{self.state} {self.score} ->\n {self.base} {self.comp}'

    def __repr__(self):
        sc = round(self.score, 3) if self.score is not None else None
        return f'<truvari.MatchResult ({self.state} {sc})>'


############################
# Parsing and set building #
############################


def file_zipper(*start_files):
    """
    Zip multiple files to yield their entries in order.

    The function takes as input tuples of (`key`, `iterable`), where:

    - `key` is an identifier (used to track which file the yielded entry comes from).
    - `iterable` is an iterable object, typically a `truvari.VariantFile`.

    The function iterates through all input files in a coordinated manner, yielding the entries in order.

    :param start_files: A variable-length argument list of tuples (`key`, `iterable`).
    :type start_files: tuple

    :yields: A tuple of (`key`, `truvari.VariantRecord`), where `key` is the file identifier and the second element is the next record from the corresponding file.
    :rtype: tuple

    :raises StopIteration: Raised when all input files have been exhausted.

    **Logs**:
        - Logs a summary of the zipping process after all files have been processed.
    """
    markers = []  # list of lists: [name, file_handler, top_entry]
    file_counts = Counter()
    for name, i in start_files:
        try:
            markers.append([name, i, next(i)])
        except StopIteration:
            # For when there are no variants in the file
            pass

    while markers:
        # Get the first entry among each of the files' tops
        first_idx = 0
        for idx, mk in enumerate(markers[1:]):
            idx += 1
            if mk[2].chrom < markers[first_idx][2].chrom or \
                    (mk[2].chrom == markers[first_idx][2].chrom and mk[2].start < markers[first_idx][2].start):
                first_idx = idx
        name, fh, entry = markers[first_idx]
        file_counts[name] += 1
        try:
            # update this file's top
            markers[first_idx][2] = next(fh)
        except StopIteration:
            # This file is done
            markers.pop(first_idx)
        yield name, entry
    logging.info("Zipped %d variants %s", sum(file_counts.values()),
                 file_counts)


def chunker(params, *files):
    """
    Given a `truvari.VariantParams` and multiple files, zip them and create chunks

    Yields tuple of the chunk of calls, and an identifier of the chunk
    """
    call_counts = Counter()
    chunk_count = 0
    cur_chrom = None
    cur_end = 0
    cur_chunk = defaultdict(list)
    reference = pysam.FastaFile(params.reference) if params.reference is not None else None
    for key, entry in file_zipper(*files):
        if entry.filter_call(key == 'base'):
            cur_chunk['__filtered'].append(entry)
            call_counts['__filtered'] += 1
            continue

        if entry.is_bnd() and params.bnddist == -1:
            cur_chunk['__filtered'].append(entry)
            call_counts['__filtered'] += 1
            continue

        if not entry.is_bnd() and entry.filter_size(key == 'base'):
            cur_chunk['__filtered'].append(entry)
            call_counts['__filtered'] += 1
            continue

        # check symbolic, resolve if needed/possible
        if entry.alts[0].startswith('<'):
            entry.resolve(reference)

        new_chrom = cur_chrom and entry.chrom != cur_chrom
        new_chunk = cur_end and cur_end + params.chunksize < entry.start
        if new_chunk or new_chrom:
            chunk_count += 1
            yield cur_chunk, chunk_count
            # Reset
            cur_chrom = None
            cur_end = 0
            cur_chunk = defaultdict(list)

        cur_chrom = entry.chrom
        cur_end = max(entry.end, cur_end)

        cur_chunk[key].append(entry)
        call_counts[key] += 1

    chunk_count += 1
    logging.info("%d chunks of %d variants %s", chunk_count,
                 sum(call_counts.values()), call_counts)
    yield cur_chunk, chunk_count
