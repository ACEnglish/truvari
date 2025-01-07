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
    A base/comp match holder
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
        # Trues are always worth more
        if self.state != other.state:
            return self.state < other.state
        m_score = self.score if self.score is not None else -float('inf')
        o_score = other.score if other.score is not None else -float('inf')
        return m_score < o_score

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
    Zip files to yield the entries in order.
    Each file must be sorted in the same order.
    start_files is a tuple of ('key', iterable)
    where key is the identifier (so we know which file the yielded entry came from)
    and iterable is usually a truvari.VariantFile

    yields key, truvari.VariantRecord
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
    unresolved_warned = False
    reference = pysam.FastaFile(params.reference) if params.reference is not None else None

    for key, entry in file_zipper(*files):
        if entry.filter_call(key == 'base'):
            cur_chunk['__filtered'].append(entry)
            call_counts['__filtered'] += 1
            continue

        if not entry.is_bnd() and entry.filter_size(key == 'base'):
            cur_chunk['__filtered'].append(entry)
            call_counts['__filtered'] += 1
            continue

        # check symbolic, resolve if needed/possible
        if params.pctseq != 0 and entry.alts[0].startswith('<'):
            was_resolved = entry.resolve(reference)
            if not was_resolved:
                if not unresolved_warned:
                    logging.warning("Some symbolic SVs couldn't be resolved")
                    unresolved_warned = True
                cur_chunk['__filtered'].append(entry)
                call_counts['__filtered'] += 1
                continue

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

        if entry.is_bnd():
            cur_chunk[f'{key}_BND'].append(entry)
        else:
            cur_chunk[key].append(entry)
        call_counts[key] += 1

    chunk_count += 1
    logging.info("%d chunks of %d variants %s", chunk_count,
                 sum(call_counts.values()), call_counts)
    yield cur_chunk, chunk_count
