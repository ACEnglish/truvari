"""
Truvari class wrapping pysam VariantRecord
"""
import re
import hashlib
import logging

import truvari
from truvari.annotations.af_calc import allele_freq_annos

RC = str.maketrans("ATCG", "TAGC")
SV_ALT_MATCH = re.compile(r"\<(?P<SVTYPE>.*)\>")

# pylint: disable=attribute-defined-outside-init


class VariantRecord:
    """
    Wrapper around pysam.VariantRecords with helper functions of variant properties and basic comparisons
    """

    def __init__(self, record, params=None):
        """
        Initialize with just the internal record
        """
        self._record = record
        self._resolved_ref = None
        self._resolved_alt = None
        self._decomp_repr = None
        self._afannos = None
        self._varsize = None
        self._vartype = None
        self.end = record.stop
        if params is None:
            self.params = truvari.VariantParams()
        else:
            self.params = params

    def __getattr__(self, name):
        """
        Delegate attribute access to the original VariantRecord
        """
        return getattr(self._record, name)

    def __setattr__(self, name, value):
        """
        Attempt to delegate attribute setting to the original VariantRecord first
        """
        if name.startswith("_") or not hasattr(self, "_record"):
            # Directly set internal attributes or during __init__ before _record is set
            super().__setattr__(name, value)
        else:
            try:
                # Try to set the attribute on the wrapped _record
                setattr(self._record, name, value)
            except AttributeError:
                # If the wrapped object does not have the attribute, set it on self
                super().__setattr__(name, value)

    def __str__(self):
        return str(self._record)

    def allele_freq_annos(self, samples=None):
        """
        Calculate allele annotations for a VCF Entry

        :param `samples`: Subset of samples from the entry over which to calculate annos
        :type `samples`: list of strings, optional

        :return: | Dictonary of
                 | AF - allele frequency
                 | MAF - minor allele frequency
                 | ExcHet - excess heterozygosity
                 | HWE - hardy weinberg equilibrium
                 | AC - allele count
                 | MAC - minor allele count
                 | AN - number of called alleles
        :rtype: dict

        Example
            >>> import truvari
            >>> v = truvari.VariantFile('repo_utils/test_files/variants/multi.vcf.gz')
            >>> e = next(v)
            >>> e.allele_freq_annos()
            {'AF': 0.5, 'MAF': 0.5, 'ExcHet': 1.0, 'HWE': 1.0, 'MAC': 1, 'AC': [1, 1], 'AN': 2, 'N_HEMI': 0, 'N_HOMREF': 0, 'N_HET': 1, 'N_HOMALT': 0, 'N_MISS': 2}
        """
        if self._afannos is None:
            self._afannos = allele_freq_annos(self, samples)
        return self._afannos

    def bnd_direction_strand(self):
        """
        Parses a BND ALT string to determine its direction and strand.

        A BND (breakend) ALT string indicates a structural variant breakpoint. This method parses the ALT string to determine:

        - The direction: "left" means the piece is anchored on the left side of the breakpoint, while "right" means it's anchored on the right.
        - The strand: "direct" indicates the base is on the direct strand, and "complement" indicates the base is on the complement strand.

        .. note::
            This method assumes that `self.is_bnd()` is `True`, meaning the variant is a BND-type structural variant.

        :return: A tuple containing the direction ("left" or "right") and the strand ("direct" or "complement").
        :rtype: tuple (str, str)

        :raises ValueError: If the ALT string does not follow the expected BND format.
        """
        bnd = self.alts[0]
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

    def bnd_position(self):
        """
        Extracts the chromosome and position from a BND ALT string.

        Breakend (BND) ALT strings indicate structural variant breakpoints that span across chromosomes or positions.
        This method parses the ALT string to extract the target chromosome and position of the breakpoint.

        :return: A tuple containing the chromosome (as a string) and the position (as an integer).
        :rtype: tuple (str, int)

        :raises ValueError: If the ALT string does not follow the expected BND format.

        """
        # Regular expression to match the BND format and extract chrom:pos
        match = re.search(r'[\[\]]([^\[\]:]+):(\d+)[\[\]]', self.alts[0])
        if not match:
            raise ValueError(f"Invalid BND ALT format: {self.alts[0]}")

        chrom = match.group(1)  # Extract the chromosome
        pos = int(match.group(2))  # Extract the position as an integer

        return chrom, pos

    def bnd_match(self, other):
        """
        Build a MatchResult for bnds
        """
        def bounds(entry, pos, key):
            """
            Inflate a bnd position based on CIPOS.
            """
            start = pos - self.params.bnddist
            end = pos + self.params.bnddist

            key = 'CI' + key
            idx = 0 if key == 'POS' else 1
            if key in entry.info:
                start -= abs(entry.info[key][idx])
                end += abs(entry.info[key][idx])

            return start, end

        ret = truvari.MatchResult()
        ret.base = self
        ret.comp = other

        # Only annotate distance if same chrom
        if self.chrom != other.chrom:
            logging.debug("%s and %s BND CHROM", str(self), str(other))
            return ret

        ret.st_dist = self.pos - other.pos
        ovl = truvari.overlaps(*bounds(self, self.pos, "POS"),
                               *bounds(other, other.pos, "POS"))
        if not ovl:
            logging.debug("%s and %s BND POS not within BNDDIST",
                          str(self), str(other))
            return ret

        b_pos2 = self.bnd_position()
        c_pos2 = other.bnd_position()
        if b_pos2[0] != c_pos2[0]:
            logging.debug("%s and %s BND join CHROM", str(self), str(other))
            return ret

        ret.ed_dist = b_pos2[1] - c_pos2[1]
        ovl = truvari.overlaps(*bounds(self, b_pos2[1], "END"),
                               *bounds(other, c_pos2[1], "END"))

        if not ovl:
            logging.debug(
                "%s and %s BND join POS not within BNDDIST", str(self), str(other))
            return ret

        b_bnd = self.bnd_direction_strand()
        c_bnd = other.bnd_direction_strand()

        ovl = b_bnd == c_bnd
        if not ovl:
            logging.debug("%s and %s BND strand/direction mismatch",
                          str(self), str(other))
            return ret

        if not self.params.skip_gt:
            self.compare_gts(other, ret)

        # Score is percent of allowed distance needed to find this match
        if self.params.bnddist > 0:
            ret.score = max(0, (1 - ((abs(ret.st_dist) + abs(ret.ed_dist)) / 2)
                                / (self.params.bnddist * 2)) * 100)
        else:
            ret.score = int(ret.state) * 100

        ret.state = True

        return ret

    def boundaries(self, ins_inflate=False):
        """
        Get the start and end of an entry

        :param `ins_inflate`: inflate INS boundaries by sv length
        :type `ins_inflate`: bool, optional

        :return: start/end
        :rtype: tuple (int, int)
        """
        start = self.start
        end = self.end
        if ins_inflate and self.var_type() == truvari.SV.INS:
            size = self.var_size()
            start -= size // 2
            end += size // 2
        return start, end

    def compare_gts(self, other, match):
        """
        Populates genotype-specific comparison details in a `MatchResult`.

        This method compares the genotypes of a "base" (reference) variant and a "comparison" variant and updates
        the provided `MatchResult` object with the results. It computes the genotype counts and determines the
        difference between the base and comparison genotypes.

        :param other: The other variant entry.
        :type other: `truvari.VariantRecord`
        :param match: The `MatchResult` object to update with genotype comparison details.
        :type match: truvari.MatchResult

        Updates the `MatchResult` object with the following attributes:
            - `base_gt`: The genotype of the base sample.
            - `base_gt_count`: The count of the reference allele (1) in the base genotype.
            - `comp_gt`: The genotype of the comparison sample.
            - `comp_gt_count`: The count of the reference allele (1) in the comparison genotype.
            - `gt_match`: The absolute difference between `base_gt_count` and `comp_gt_count`.
        """
        b_gt = self.gt(self.params.bSample)
        c_gt = other.gt(self.params.cSample)
        if b_gt:
            match.base_gt = b_gt
            match.base_gt_count = sum(1 for _ in match.base_gt if _ == 1)
        if c_gt:
            match.comp_gt = c_gt
            match.comp_gt_count = sum(1 for _ in match.comp_gt if _ == 1)
        match.gt_match = abs(match.base_gt_count - match.comp_gt_count)

    def copy(self):
        """
        Copies this truvari.VariantRecord
        """
        return VariantRecord(self._record.copy(), params=self.params)

    def cpx_match(self, other):
        """
        Attempt to match non-BND variants with a BND.
        This is performed by decomposing the non-BND into its BND representations before calling bnd_match.
        This will make at least two BNDs, so the MatchResults are sorted and the greater is returned.
        """
        if self.is_bnd():
            decomp = other.decompose()
            matches = [self.bnd_match(d) for d in decomp]
        else:
            decomp = self.decompose()
            matches = [d.bnd_match(other) for d in decomp]

        # Not all variants can be decomposed
        if not matches:
            mat = truvari.MatchResult()
        else:
            mat = sorted(matches)[-1]
        mat.base = self
        mat.comp = other
        return mat

    def decompose(self):
        """
        Decompose SVs into BND Records

        This method decomposes  structural variants (SVs) of `<DEL>`, `<DUP>`, and `<INV>` into
        Breakend (BND) records. The decomposed variants are stored internally and returned as a list.

        Returns:
            list: A list of new BND variant records. List will be empty when no valid decomposition

        Notes:
            - For `INV`, four BNDs are created.
            - For `DEL` and `DUP`, two BNDs are created.
            - `DUP` are assumed to be of type `DUP:TANDEM`.
        """
        # No need to make twice
        if self._decomp_repr is not None:
            return self._decomp_repr

        svtype = self.var_type()
        ret = []

        if svtype == truvari.SV.INV:
            record1 = self.copy()
            record1.alts = (f"N]{self.chrom}:{self.end}]",)
            record1.info["SVTYPE"] = "BND"

            record2 = self.copy()
            record2.pos += 1
            record2.alts = (f"[{self.chrom}:{self.end + 1}[N",)
            record2.info["SVTYPE"] = "BND"

            record3 = self.copy()
            record3.pos = self.end
            record3.alts = (f"N]{self.chrom}:{self.pos}]",)
            record3.info["SVTYPE"] = "BND"

            record4 = self.copy()
            record4.pos = self.end + 1
            record4.alts = (f"[{self.chrom}:{self.pos + 1}[N",)
            record4.info["SVTYPE"] = "BND"

            ret = [record1, record2, record3, record4]

        elif svtype == truvari.SV.DEL:
            record1 = self.copy()
            record1.alts = (f"N[{self.chrom}:{self.end + 1}[",)
            record1.info["SVTYPE"] = "BND"

            record2 = self.copy()
            record2.pos = self.end + 1
            record2.alts = (f"]{self.chrom}:{self.pos}]N",)
            record2.info["SVTYPE"] = "BND"

            ret = [record1, record2]

        elif svtype == truvari.SV.DUP:
            # Assumes DUP:TANDEM
            record1 = self.copy()
            record1.pos += 1
            record1.alts = (f"]{self.chrom}:{self.end}]N",)
            record1.info["SVTYPE"] = "BND"

            record2 = self.copy()
            record2.pos = self.end
            record2.alts = (f"N[{self.chrom}:{self.pos + 1}[",)
            record2.info["SVTYPE"] = "BND"

            ret = [record1, record2]

        self._decomp_repr = ret
        return ret

    def distance(self, other):
        """
        Calculate the start and end distances of the pair. Negative distances
        indicate self is upstream of other

        :param `other`: Other to compare
        :type `other`: :class:`truvari.VariantRecord`

        :return: starts and ends distance
        :rtype: tuple (int, int)
        """
        astart, aend = self.boundaries()
        bstart, bend = other.boundaries()
        return astart - bstart, aend - bend

    def gt(self, sample=0):
        """
        Returns the raw GT field from a record for a specific sample.
        Returns None if GT is not in FORMAT
        """
        if "GT" not in self.samples[sample]:
            return None
        return self.samples[sample]['GT']

    def get_ref(self):
        """
        Returns either the resolved REF or _record.ref
        """
        if self._resolved_ref:
            return self._resolved_ref
        return self._record.ref

    def get_alt(self):
        """
        Returns either the resolved ALT or _record.alt[0]
        """
        if self._resolved_alt:
            return self._resolved_alt
        return self._record.alts[0]

    def is_bnd(self):
        """
        Returns if a record is a resolved BND
        """
        return '[' in self._record.alts[0] or ']' in self._record.alts[0]

    def is_filtered(self, values=None):
        """
        Checks if entry should be filtered given the filter values provided.
        If `values` is None, assume that filter must have PASS or be blank '.'

        :param `values`: set of filter values for intersection
        :type `values`: set, optional

        :return: True if entry should be filtered
        :rtype: bool

        Example
            >>> import truvari
            >>> v = truvari.VariantFile('repo_utils/test_files/variants/input1.vcf.gz')
            >>> e = next(v)
            >>> e.is_filtered() # PASS shouldn't be filtered
            False
            >>> e = next(v)
            >>> e.is_filtered(["lowQ"]) # Call isn't lowQ, so filter
            True
        """
        if values is None:
            return len(self.filter) != 0 and 'PASS' not in self.filter
        return len(set(values).intersection(set(self.filter))) == 0

    def is_monrefstar(self):
        """
        Returns if a record is a reference monozygotic or star-alt
        """
        return not self.alts or self.alts[0] in (None, '*')

    def is_multi(self):
        """
        Returns if a record is multi-allelic
        """
        return len(self.alts) > 1

    def is_present(self, sample=0, allow_missing=True):
        """
        Checks if entry's sample genotype is present and is heterozygous or homozygous (a.k.a. present)
        If allow_missing, just check for a 1 in the genotype. Otherwise, a missing ('.') genotype isn't
        considered present

        :param `sample`: sample name
        :type `sample`: string, optional

        :return: True if variant is present in the sample
        :rtype: bool

        Example
            >>> import truvari
            >>> v = truvari.VariantFile('repo_utils/test_files/variants/input1.vcf.gz')
            >>> e = next(v)
            >>> e.is_present(allow_missing=False)
            True
        """
        gt = self.gt(sample)
        if allow_missing:
            return 1 in gt
        return truvari.get_gt(gt) in [truvari.GT.HET, truvari.GT.HOM]

    def is_single_bnd(self):
        """
        Returns if a record is a single-end BND
        """
        return '.' in (self.alts[0][0], self.alts[0][-1])

    def is_symbolic(self):
        """
        Returns if a record is a symbolic variant (e.g. ALT == <DEL>)
        """
        return self.alts[0][0] == '<'

    def is_resolved(self):
        """
        Checks if an SV is a sequence resolved by checking it is
        has an alt, is not symbolic, is not monref, is not BND

        :return: is_resolved
        :rtype: bool
        """
        orig_resolved = not self.is_monrefstar() and \
            not self.is_symbolic() and \
            not self.is_bnd()
        truv_resolved = self._resolved_ref is not None \
            and self._resolved_alt is not None
        return orig_resolved or truv_resolved

    def filter_call(self, base=False):
        """
        Determines whether a variant call should be filtered based on Truvari parameters or specific requirements.

        This method evaluates a variant entry (`entry`) and checks if it should be excluded from further processing
        based on filtering criteria such as monomorphic reference, multi-allelic records, filtering status,
        sample presence, or unsupported single-end BNDs.

        :param entry: The variant entry to evaluate.
        :type entry: truvari.VariantRecord
        :param base: A flag indicating whether the entry is the "base" (reference) call or the "comparison" call.
                     Filtering behavior may differ based on this flag.
        :type base: bool, optional

        :return: `True` if the variant should be filtered (excluded), otherwise `False`.
        :rtype: bool

        :raises ValueError: If the entry is multi-allelic and `check_multi` is enabled in the Truvari parameters.

        Filtering Logic:
            - **Monomorphic Reference:** If `check_monref` is enabled and the entry is a monomorphic reference, it is filtered.
            - **Multi-Allelic Records:** If `check_multi` is enabled and the entry is multi-allelic, an error is raised.
            - **Filtered Variants:** If `passonly` is enabled and the entry is flagged as filtered, it is excluded.
            - **Sample Presence:** If `no_ref` is set to include the entry's type (base or comparison) or `pick == 'ac'`, the sample must be present in the entry.
            - **Single-End BNDs:** Single-end BNDs are always excluded.
        """
        if self.params.check_monref and self.is_monrefstar():
            return True

        if self.params.check_multi and self.is_multi():
            raise ValueError(
                f"Cannot compare multi-allelic records. Please split\nline {str(self)}")

        if self.params.passonly and self.is_filtered():
            return True

        prefix = 'b' if base else 'c'
        if self.params.no_ref in ["a", prefix] or self.params.pick == 'ac':
            samp = self.params.bSample if base else self.params.cSample
            if not self.is_present(samp):
                return True

        if self.params.no_single_bnd and self.is_single_bnd():
            return True
        return False

    def filter_size(self, base=False):
        """
        Determines whether a variant entry should be filtered based on its size.

        This method evaluates the size of a variant and checks if it falls outside the specified size thresholds.
        Filtering criteria depend on whether the entry is a "base" (reference) or "comparison" call.

        :param entry: The variant entry to evaluate.
        :type entry: truvari.VariantRecord
        :param base: A flag indicating whether the entry is the "base" (reference) call. If `True`, the `sizemin`
                     parameter is used as the minimum size threshold. Otherwise, the `sizefilt` parameter is used.
        :type base: bool, optional

        :return: `True` if the variant should be filtered due to its size, otherwise `False`.
        :rtype: bool

        Filtering Logic:
            - **Minimum Size (Base):** If `base=True` and the size is less than `sizemin`, the variant is filtered.
            - **Minimum Size (Comparison):** If `base=False` and the size is less than `sizefilt`, the variant is filtered.
            - **Maximum Size:** If the size exceeds `sizemax`, the variant is filtered.
        """
        size = self.var_size()
        return (self.params.sizemax != -1 and size > self.params.sizemax) \
            or (base and size < self.params.sizemin) \
            or (not base and size < self.params.sizefilt)

    def get_record(self):
        """
        Return internal pysam.VariantRecord
        """
        return self._record

    def match(self, other):
        """
        Build a MatchResult from comparison of two VariantRecords
        If self and other are non-bnd, calls VariantRecord.var_match,
        If self and other are bnd, calls VariantRecord.bnd_match
        If decompose is on, calls VariantRecord.cpx_match
        Otherwise, returns a default MatchResult
        """
        if not self.is_bnd() and not other.is_bnd():
            return self.var_match(other)
        if self.is_bnd() and other.is_bnd():
            return self.bnd_match(other)
        if self.params.decompose:
            return self.cpx_match(other)
        mat = truvari.MatchResult()
        mat.base = self
        mat.comp = other
        return mat

    def move_record(self, out_vcf, take_samples=None, into_samples=None):
        """
        Similar to pysam.VariantRecord.translate, except this allows adding new
        FORMAT fields and optionally selects a subset of samples
        """
        ret = out_vcf.new_record(contig=self.contig, start=self.start, stop=self.stop, alleles=self.alleles,
                                 id=self.id, qual=self.qual, filter=self.filter, info=self.info)
        if take_samples is None:
            take_samples = self.samples.keys()
        if into_samples is None:
            into_samples = list(range(len(take_samples)))
        for t, i in zip(take_samples, into_samples):
            for k, v in self.samples[t].items():
                try:
                    ret.samples[i][k] = v
                except TypeError:
                    # Somtimes they're supposed to be tuples.
                    pass
        return VariantRecord(ret)

    def recovl(self, other, ins_inflate=True):
        """
        Calculates reciprocal overlap of two entries

        :param `other`: Other entry to compare with
        :type `other`: :class:`truvari.VariantRecord`
        :param `ins_inflate`: inflate entry boundaries
        :type `ins_inflate`: bool

        :return: The reciprocal overlap
        :rtype: float

        Example
            >>> import truvari
            >>> v = truvari.VariantFile('repo_utils/test_files/variants/input1.vcf.gz')
            >>> a = next(v)
            >>> b = next(v)
            >>> a.recovl(b)
            0
        """
        astart, aend = self.boundaries(ins_inflate)
        bstart, bend = other.boundaries(ins_inflate)
        return truvari.reciprocal_overlap(astart, aend, bstart, bend)

    def resolve(self, ref):
        """
        Attempts to resolve the REF/ALT sequences of a structural variant (SV) and updates the instance attributes.

        This method tries to reconstruct the REF and ALT sequences for specific types of structural variants
        (e.g., deletions, inversions, duplications) based on the provided reference genome sequence. If successful,
        it stores the resolved sequences in `self._resolved_ref` and `self._resolved_alt`, and may adjust `self.end`.

        :param ref: The reference genome used to resolve the structural variant.
                    Typically a `pysam.FastaFile` or equivalent object providing access to reference sequences.
        :type ref: pysam.FastaFile

        :return: `True` if the REF/ALT sequences are successfully resolved, otherwise `False`.
        :rtype: bool

        Updates:
            - `self._resolved_ref`: The resolved REF sequence.
            - `self._resolved_alt`: The resolved ALT sequence.
            - `self.end`: May be adjusted for specific variant types (e.g., duplications).

        Resolution Logic:
            - If `ref` is `None`, the ALT is `<CNV>` or `<INS>`, or the start position is invalid, the method returns `False`.
            - The method fetches the sequence from the reference genome using the variant's `chrom`, `start`, and `end`.
            - The resolution depends on the variant type:
                - **Deletion (`DEL`)**: The REF sequence is the fetched reference sequence, and the ALT is the first base.
                - **Inversion (`INV`)**: The REF sequence is the fetched reference sequence, and the ALT is the reverse complement of the REF.
                - **Duplication (`DUP`)** (when `self.params.dup_to_ins` is enabled): The REF is the first base of the fetched sequence, and the ALT is the full sequence. The `end` is adjusted to be `start + 1`.
            - Returns `False` if the variant type is unsupported or cannot be resolved.
        """
        # Can't resolve everything
        if ref is None or self.is_monrefstar() \
                or not self.is_symbolic() \
                or self.start > ref.get_reference_length(self.chrom) \
                or self.var_size() > self.params.max_resolve:
            return False

        seq = ref.fetch(self.chrom, self.start, self.end)
        svtype = self.var_type()
        if svtype == truvari.SV.DEL:
            self._resolved_ref = seq
            self._resolved_alt = seq[0]
        elif svtype == truvari.SV.INV:
            self._resolved_ref = seq
            self._resolved_alt = seq.translate(RC)[::-1]
        elif svtype == truvari.SV.DUP and self.params.dup_to_ins:
            self._resolved_ref = seq[0]
            self._resolved_alt = seq
            self.end = self.start + 1
        else:
            return False

        return True

    def same_type(self, other):
        """
        Check if self.var_type() == other.var_type() with extra handling for dup-to-ins

        :param `other`: Other entry to compare with
        :type `other`: :class:`truvari.VariantRecord`
        :param `dup_to_ins`: Convert DUP to INS types
        :type `dup_to_ins`: bool

        :return: True if entry SVTYPEs match
        :rtype: bool
        """
        a_type = self.var_type()
        b_type = other.var_type()
        if self.params.dup_to_ins and a_type == truvari.SV.DUP:
            a_type = truvari.SV.INS
        if self.params.dup_to_ins and b_type == truvari.SV.DUP:
            b_type = truvari.SV.INS
        return a_type == b_type

    def seqsim(self, other):
        """
        Calculate sequence similarity of two entries. If reference is not None,
        compare their shared reference context. Otherwise, use the unroll technique.

        Assumes that is either a sequence resolved variant or self.resolve has been called.

        :param `other`: Other to compare with
        :type `other`: :class:`truvari.VariantRecord`
        :param `roll`: compare the lexicographically minimum rotation of sequences
        :type `roll`: bool

        :return: sequence similarity
        :rtype: float
        """
        # Shortcut to save compute - probably unneeded
        if self.get_ref() == other.get_ref() and self.get_alt() == other.get_alt():
            return 1.0

        # Inversions aren't rolled
        if (self.var_type() == truvari.SV.INV and other.var_type() == truvari.SV.INV):
            allele1 = self.get_alt()
            allele2 = other.get_alt()
            return truvari.seqsim(allele1, allele2)

        a_seq = self.get_ref() if self.var_type() == truvari.SV.DEL else self.get_alt()
        a_seq = a_seq.upper()
        b_seq = other.get_ref() if other.var_type() == truvari.SV.DEL else other.get_alt()
        b_seq = b_seq.upper()
        st_dist, ed_dist = self.distance(other)

        if self.params.no_roll or st_dist == 0 or ed_dist == 0:
            return truvari.seqsim(a_seq, b_seq)

        # Return best of rolled, unrolled from both ends, and direct similarity
        # Whichever is highest is how similar these sequences can be
        return truvari.best_seqsim(a_seq, b_seq, st_dist)

    def sizesim(self, other):
        """
        Calculate the size similarity and difference for two entries

        :param `other`: Other entry to compare with
        :type `other`: :class:`truvari.VariantRecord`

        :return: size similarity and size diff (A - B)
        :rtype: (float, int)

        Example
            >>> import truvari
            >>> v = truvari.VariantFile('repo_utils/test_files/variants/input1.vcf.gz')
            >>> a = next(v)
            >>> b = next(v)
            >>> a.sizesim(b)
            (0.07142857142857142, 13)
        """
        return truvari.sizesim(self.var_size(), other.var_size())

    def var_match(self, other):
        """
        Build a MatchResult
        """
        ret = truvari.MatchResult()
        ret.base = self
        ret.comp = other
        ret.state = True

        if not self.params.typeignore and not self.same_type(other):
            logging.debug("%s and %s are not the same SVTYPE",
                          str(self), str(other))
            ret.state = False
            if self.params.short_circuit:
                return ret

        bstart, bend = self.boundaries()
        cstart, cend = other.boundaries()
        if not truvari.overlaps(bstart - self.params.refdist, bend + self.params.refdist, cstart, cend):
            logging.debug("%s and %s are not within REFDIST",
                          str(self), str(other))
            ret.state = False
            if self.params.short_circuit:
                return ret

        ret.sizesim, ret.sizediff = self.sizesim(other)
        if ret.sizesim < self.params.pctsize:
            logging.debug("%s and %s size similarity is too low (%.3f)",
                          str(self), str(other), ret.sizesim)
            ret.state = False
            if self.params.short_circuit:
                return ret

        if not self.params.skip_gt:
            self.compare_gts(other, ret)

        ret.ovlpct = self.recovl(other)
        if ret.ovlpct < self.params.pctovl:
            logging.debug("%s and %s overlap percent is too low (%.3f)",
                          str(self), str(other), ret.ovlpct)
            ret.state = False
            if self.params.short_circuit:
                return ret

        if self.params.pctseq > 0 and self.is_resolved() and other.is_resolved():
            ret.seqsim = self.seqsim(other)
            if ret.seqsim < self.params.pctseq:
                logging.debug("%s and %s sequence similarity is too low (%.3ff)",
                              str(self), str(other), ret.seqsim)
                ret.state = False
                if self.params.short_circuit:
                    return ret
        else:
            ret.seqsim = 0

        ret.st_dist, ret.ed_dist = bstart - cstart, bend - cend
        ret.calc_score()

        return ret

    def var_size(self):
        """
        Determine the size of the variant.

        .. note:: How size is determined

            - Starts by trying to use INFO/SVLEN
            - If SVLEN is unavailable and ALT field is an SV (e.g. <INS>, <DEL>, etc), \
            use abs(vcf.start - vcf.end). The INFO/END tag needs to be available, \
            especially for INS.
            - If len of vcf.REF and vcf.ALT[0] are equal, return len of vcf.REF
            - Otherwise, return the size difference of the sequence resolved call using \
            abs(len(vcf.REF) - len(str(vcf.ALT[0])))

        :return: the entry's size
        :rtype: int
        """
        # Visited enough that we take the memory hit
        if self._varsize is not None:
            return self._varsize

        if "SVLEN" in self.info:
            if type(self.info["SVLEN"]) in [list, tuple]:
                size = self.info["SVLEN"][0]
            else:
                size = self.info["SVLEN"]
            try:
                size = abs(int(size))
                self._varsize = size
                return size
            except ValueError:
                pass

        if not self.is_resolved() and self.alts is not None and self.alts[0].count("<"):
            start, end = self.boundaries()
            size = end - start
        else:
            r_len = len(self.get_ref())
            a_len = len(self.get_alt())
            if r_len == a_len:
                if r_len == 1:
                    size = 0  # SNPs are special
                else:
                    size = r_len
            else:
                size = abs(r_len - a_len)
        self._varsize = size
        return size

    def var_type(self):
        """
        Return entry's SVTYPE

        .. note::
            the first of these checks that succeeds determins the svtype

            - BND if `pysam.VariantRecord.alleles_variant_types[1]` == "BND"
            - Checking INFO/SVTYPE
            - if `is_resolved`, REF and ALT seq size difference
            - if symbolic (e.g. <DEL>) parse the ALT. Subtypes are ignored \
            (e.g. `DUP:ME` becomes `DUP`)
            - svtype is unknown (UNK)

        :return: SV type
        :rtype: :class:`truvari.SV`
        """
        if self._vartype is not None:
            return self._vartype
        try:
            if self.alleles_variant_types[1] == 'BND':
                self._vartype = truvari.SV.BND
                return self._vartype
        except (IndexError, AttributeError):
            pass

        if "SVTYPE" in self.info:
            ret_type = self.info["SVTYPE"]
            if isinstance(ret_type, (list, tuple)):
                ret_type = ret_type[0]
            
            self._vartype = truvari.get_svtype(ret_type)
            return self._vartype

        if self.is_resolved():
            if len(self.ref) < len(self.alts[0]):
                ret_type = "INS"
            elif len(self.ref) > len(self.alts[0]):
                ret_type = "DEL"
            elif len(self.ref) == len(self.alts[0]):
                ret_type = "SNP" if len(self.ref) == 1 else "UNK"
            self._vartype = truvari.get_svtype(ret_type)
            return self._vartype

        mat = SV_ALT_MATCH.match(
            self.alts[0]) if self.alts is not None else None
        if mat is not None:
            self._vartype = truvari.get_svtype(mat.groupdict()["SVTYPE"].split(':')[0])
            return self._vartype

        self._vartype = truvari.get_svtype("UNK")
        return self._vartype

    def to_hash(self, hasher=hashlib.sha1):
        """
        Turn variant into a key and hash with provided hasher

        :param `hasher`: hashing function
        :type `entry`: method

        :return: hash
        :rtype: string
        """
        return hasher(self.to_key().encode()).hexdigest()

    def to_key(self, prefix="", bounds=False):
        """
        Turn a vcf entry into a hashable key string. Use the prefix (base/comp) to indicate the source
        VCF when consolidating multiple files' calls. If bounds: call entry_boundaries for start/stop.

        .. warning::
            If a caller redundantly calls a variant exactly the same it will not have a unique key

        :param `prefix`: prefix
        :type `prefix`: string, optional
        :param `bounds`: use entry_boundaries
        :type `bounds`: bool, optional

        :return: hashable string uniquely identifying the variant
        :rtype: string
        """
        if prefix:
            prefix += '.'
        alt = self.alts[0] if self.alts is not None else "."
        if bounds:
            start, end = self.boundaries()
            return f"{prefix}{self.chrom}:{start}-{end}({self.ref}|{alt})"
        return f"{prefix}{self.chrom}:{self.start}-{self.end}.{alt}"

    def within_tree(self, tree):
        """
        Extract entry and tree boundaries to call `truvari.coords_within`
        """
        qstart, qend = self.boundaries()
        if self.chrom not in tree:
            return False
        m_ovl = tree[self.chrom].overlap(qstart, qend)
        if len(m_ovl) != 1:
            return False
        m_ovl = list(m_ovl)[0]
        end_within = self.var_type() != truvari.SV.INS

        return truvari.coords_within(qstart, qend, m_ovl.begin, m_ovl.end - 1, end_within)

    def within(self, rstart, rend):
        """
        Extract entry boundaries and type to call `truvari.coords_within`
        """
        qstart, qend = self.boundaries()
        end_within = self.var_type() != truvari.SV.INS
        return truvari.coords_within(qstart, qend, rstart, rend, end_within)
