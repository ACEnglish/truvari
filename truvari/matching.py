"""
Comparison engine
"""
import re
import types
import logging
from collections import Counter, defaultdict
from functools import total_ordering
import pysam
import truvari


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


class Matcher():
    """
    Holds matching parameters. Allows calls to be checked for filtering and matches to be made

    Example
        >>> import pysam
        >>> import truvari
        >>> mat = truvari.Matcher()
        >>> mat.params.pctseq = 0
        >>> v = pysam.VariantFile('repo_utils/test_files/variants/input1.vcf.gz')
        >>> one = next(v); two = next(v)
        >>> mat.build_match(one, two)
        <truvari.bench.MatchResult (False 2.381)>

    Look at `Matcher.make_match_params()` for a list of all params and their defaults
    """

    def __init__(self, args=None):
        """
        Initalize. args is a Namespace from argparse
        """
        if args is not None:
            self.params = self.make_match_params_from_args(args)
        else:
            self.params = self.make_match_params()

        self.reference = None
        if self.params.reference is not None:
            self.reference = pysam.FastaFile(self.params.reference)

    @staticmethod
    def make_match_params():
        """
        Makes a simple namespace of matching parameters. Holds defaults
        """
        params = types.SimpleNamespace()
        params.reference = None
        params.refdist = 500
        params.pctseq = 0.70
        params.pctsize = 0.70
        params.pctovl = 0.0
        params.typeignore = False
        params.no_roll = False
        params.chunksize = 1000
        params.bSample = 0
        params.cSample = 0
        params.dup_to_ins = False
        params.bnddist = 100
        params.sizemin = 50
        params.sizefilt = 30
        params.sizemax = 50000
        params.passonly = False
        params.no_ref = False
        params.pick = 'single'
        params.ignore_monref = True
        params.check_multi = True
        params.check_monref = True
        return params

    @staticmethod
    def make_match_params_from_args(args):
        """
        Makes a simple namespace of matching parameters
        """
        ret = types.SimpleNamespace()
        ret.reference = args.reference
        ret.refdist = args.refdist
        ret.pctseq = args.pctseq
        ret.pctsize = args.pctsize
        ret.pctovl = args.pctovl
        ret.typeignore = args.typeignore
        ret.no_roll = args.no_roll
        ret.chunksize = args.chunksize
        ret.bSample = args.bSample if args.bSample else 0
        ret.cSample = args.cSample if args.cSample else 0
        ret.dup_to_ins = args.dup_to_ins if "dup_to_ins" in args else False
        ret.bnddist = args.bnddist if 'bnddist' in args else -1
        # filtering properties
        ret.sizemin = args.sizemin
        ret.sizefilt = args.sizefilt
        ret.sizemax = args.sizemax
        ret.passonly = args.passonly
        ret.no_ref = args.no_ref
        ret.pick = args.pick if "pick" in args else "single"
        ret.check_monref = True
        ret.check_multi = True
        return ret

    def filter_call(self, entry, base=False):
        """
        Returns True if the call should be filtered based on parameters or truvari requirements
        Base has different filtering requirements, so let the method know
        """
        if self.params.check_monref and entry.is_monrefstar():
            return True

        if self.params.check_multi and entry.is_multi():
            raise ValueError(
                f"Cannot compare multi-allelic records. Please split\nline {str(entry)}")

        if self.params.passonly and entry.is_filtered():
            return True

        prefix = 'b' if base else 'c'
        if self.params.no_ref in ["a", prefix] or self.params.pick == 'ac':
            samp = self.params.bSample if base else self.params.cSample
            if not entry.is_present(samp):
                return True

        # No single end BNDs
        return entry.is_single_bnd()

    def size_filter(self, entry, base=False):
        """
        Returns True if entry should be filtered due to its size
        """
        size = entry.size()
        return (size > self.params.sizemax) \
            or (base and size < self.params.sizemin) \
            or (not base and size < self.params.sizefilt)

    def compare_gts(self, match, base, comp):
        """
        Given a MatchResult, populate the genotype specific comparisons in place
        """
        b_gt = base.gt(self.params.bSample)
        c_gt = comp.gt(self.params.cSample)
        if b_gt:
            match.base_gt = b_gt
            match.base_gt_count = sum(1 for _ in match.base_gt if _ == 1)
        if c_gt:
            match.comp_gt = c_gt
            match.comp_gt_count = sum(1 for _ in match.comp_gt if _ == 1)
        match.gt_match = abs(match.base_gt_count - match.comp_gt_count)

    def build_match(self, base, comp, matid=None, skip_gt=False, short_circuit=False):
        """
        Build a MatchResult
        if skip_gt, don't do genotype comparison
        if short_circuit, return after first failure
        """
        ret = MatchResult()
        ret.base = base
        ret.comp = comp

        ret.matid = matid
        ret.state = True

        if not self.params.typeignore and not base.same_type(comp, self.params.dup_to_ins):
            logging.debug("%s and %s are not the same SVTYPE",
                          str(base), str(comp))
            ret.state = False
            if short_circuit:
                return ret

        bstart, bend = base.boundaries()
        cstart, cend = comp.boundaries()
        if not truvari.overlaps(bstart - self.params.refdist, bend + self.params.refdist, cstart, cend):
            logging.debug("%s and %s are not within REFDIST",
                          str(base), str(comp))
            ret.state = False
            if short_circuit:
                return ret

        ret.sizesim, ret.sizediff = base.sizesim(comp)
        if ret.sizesim < self.params.pctsize:
            logging.debug("%s and %s size similarity is too low (%.3f)",
                          str(base), str(comp), ret.sizesim)
            ret.state = False
            if short_circuit:
                return ret

        if not skip_gt:
            self.compare_gts(ret, base, comp)

        ret.ovlpct = base.recovl(comp)
        if ret.ovlpct < self.params.pctovl:
            logging.debug("%s and %s overlap percent is too low (%.3f)",
                          str(base), str(comp), ret.ovlpct)
            ret.state = False
            if short_circuit:
                return ret

        if self.params.pctseq > 0:
            ret.seqsim = base.seqsim(comp, self.params.no_roll)
            if ret.seqsim < self.params.pctseq:
                logging.debug("%s and %s sequence similarity is too low (%.3ff)",
                              str(base), str(comp), ret.seqsim)
                ret.state = False
                if short_circuit:
                    return ret
        else:
            ret.seqsim = 0

        ret.st_dist, ret.ed_dist = bstart - cstart, bend - cend
        ret.calc_score()

        return ret

    def bnd_build_match(self, base, comp, matid=None, **_kwargs):
        """
        Build a MatchResult for bnds
        """
        def bounds(pos):
            """
            Inflate a bnd position based on CIPOS.
            Experimental, but it's getting tricky, so need clarification
            For example, why does GIAB use CIPOS1 instead of the standard CIPOS? (not to mention Type=String).
            Also, I assumed that the CIPOS was ±POS, but it seems like the standard is that the ambiguity is between
            POS, POS+CIPOS[0]. But then, there is no enforcement of strands with CIPOS, I think, so it could always be
            positive and if it's a complement BND it might need POS - CIPOS?
            I'm dropping for now
            entry and key were also parameters
            key = 'CI' + key
            idx = 0 if key == 'POS' else 1
            if key in entry.info: # Just CIPOS
                start -= entry.info[key][idx]
                end += entry.info[key][idx]
            else: # Special CIPOS1/CIPOS2
                if key + '1' in entry.info:
                    k = entry.info[key + '1']
                    if k not in [None, '.']:
                        start -= int(k)
                if key + '2' in entry.info:
                    k = entry.info[key + '2']
                    if k not in [None, '.']:
                        end += int(k)
            """
            start = pos - self.params.bnddist
            end = pos + self.params.bnddist

            return start, end

        ret = truvari.MatchResult()
        ret.base = base
        ret.comp = comp

        ret.matid = matid
        # Only put start distance same chrom pos2
        if base.chrom != comp.chrom:
            logging.debug("%s and %s BND CHROM", str(base), str(comp))
            return ret

        ret.st_dist = base.pos - comp.pos
        ovl = truvari.overlaps(*bounds(base.pos), *bounds(comp.pos))
        if not ovl:
            logging.debug("%s and %s BND POS not within BNDDIST",
                          str(base), str(comp))
            return ret

        b_pos2 = bnd_position(base.alts[0])
        c_pos2 = bnd_position(comp.alts[0])
        if b_pos2[0] != c_pos2[0]:
            logging.debug("%s and %s BND join CHROM", str(base), str(comp))
            return ret

        ret.ed_dist = b_pos2[1] - c_pos2[1]
        ovl = truvari.overlaps(*bounds(b_pos2[1]), *bounds(c_pos2[1]))

        if not ovl:
            logging.debug(
                "%s and %s BND join POS not within BNDDIST", str(base), str(comp))
            return ret

        b_bnd = bnd_direction_strand(base.alts[0])
        c_bnd = bnd_direction_strand(comp.alts[0])

        ovl = b_bnd == c_bnd
        if not ovl:
            logging.debug("%s and %s BND strand/direction mismatch",
                          str(base), str(comp))
            return ret

        self.compare_gts(ret, base, comp)

        # Score is percent of allowed distance needed to find this match
        if self.params.bnddist > 0:
            ret.score = max(0, (1 - ((abs(ret.st_dist) + abs(ret.ed_dist)) / 2)
                                / self.params.bnddist) * 100)
        else:
            ret.score = int(ret.state) * 100

        ret.state = True

        return ret

############################
# Parsing and set building #
############################


def file_zipper(*start_files):
    """
    Zip files to yield the entries in order.
    Each file must be sorted in the same order.
    start_files is a tuple of ('key', iterable)
    where key is the identifier (so we know which file the yielded entry came from)
    and iterable is usually a pysam.VariantFile

    yields key, pysam.VariantRecord
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
        yield name, truvari.VariantRecord(entry)
    logging.info("Zipped %d variants %s", sum(
        file_counts.values()), file_counts)


def chunker(matcher, *files):
    """
    Given a Matcher and multiple files, zip them and create chunks

    Yields tuple of the chunk of calls, and an identifier of the chunk
    """
    call_counts = Counter()
    chunk_count = 0
    cur_chrom = None
    cur_end = 0
    cur_chunk = defaultdict(list)
    unresolved_warned = False
    for key, entry in file_zipper(*files):
        if matcher.filter_call(entry, key == 'base'):
            cur_chunk['__filtered'].append(entry)
            call_counts['__filtered'] += 1
            continue

        if not entry.is_bnd() and matcher.size_filter(entry, key == 'base'):
            cur_chunk['__filtered'].append(entry)
            call_counts['__filtered'] += 1
            continue

        # check symbolic, resolve if needed/possible
        if matcher.params.pctseq != 0 and entry.alts[0].startswith('<'):
            was_resolved = entry.resolve(matcher.reference, matcher.params.dup_to_ins)
            if not was_resolved:
                if not unresolved_warned:
                    logging.warning("Some symbolic SVs couldn't be resolved")
                    unresolved_warned = True
                cur_chunk['__filtered'].append(entry)
                call_counts['__filtered'] += 1
                continue

        new_chrom = cur_chrom and entry.chrom != cur_chrom
        new_chunk = cur_end and cur_end + matcher.params.chunksize < entry.start
        if new_chunk or new_chrom:
            chunk_count += 1
            yield cur_chunk, chunk_count
            # Reset
            cur_chrom = None
            cur_end = 0
            cur_chunk = defaultdict(list)

        cur_chrom = entry.chrom
        cur_end = max(entry.stop, cur_end)

        if entry.is_bnd():
            cur_chunk[f'{key}_BND'].append(entry)
        else:
            cur_chunk[key].append(entry)
        call_counts[key] += 1

    chunk_count += 1
    logging.info("%d chunks of %d variants %s", chunk_count,
                 sum(call_counts.values()), call_counts)
    yield cur_chunk, chunk_count

####################
# Helper functions #
####################


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
