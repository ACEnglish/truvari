"""
Truvari functions wrapping pysam VariantFile and VariantRecords
"""
import re
import hashlib
import logging
import pysam
import truvari

RC = str.maketrans("ATCG", "TAGC")


class VariantFile:
    """
    Wrapper around pysam.VariantFile with helper functions for iteration
    Note: The context manager functionality of pysam.VariantFile is not available with truvari.VariantFile
    """

    def __init__(self, filename, *args, **kwargs):
        self._vcf = pysam.VariantFile(filename, *args, **kwargs)

    def __getattr__(self, name):
        """
        Delegate attribute access to the original VariantFile
        """
        return getattr(self._vcf, name)

    def __iter__(self):
        """
        Iterate the VariantFile, wrapping into truvari VariantRecords
        """
        for i in self._vcf:
            yield VariantRecord(i)

    def __next__(self):
        """
        Return the next
        """
        return VariantRecord(next(self._vcf))

    def fetch(self, *args, **kwargs):
        """
        Fetch from the VariantFile, wrapping into truvari VariantRecords
        """
        for i in self._vcf.fetch(*args, **kwargs):
            yield truvari.VariantRecord(i)

    def write(self, record, resolved=False):
        """
        Pull pysam VarianRecord out of truvari VariantRecord before writing
        If resolved, replace the record's REF and ALT with self.get_ref() self.get_alt()
        """
        out = record.get_record()
        if resolved:
            out.ref = record.get_ref()
            out.alts = (record.get_alt(),)
        self._vcf.write(out)


class VariantRecord:
    """
    Wrapper around pysam.VariantRecords with helper functions of variant properties and basic comparisons
    """

    def __init__(self, record):
        """
        Initialize with just the internal record
        """
        self._record = record
        self.resolved_ref = None
        self.resolved_alt = None
        self.end = record.stop

    def __getattr__(self, name):
        """
        Delegate attribute access to the original VariantRecord
        """
        return getattr(self._record, name)

    def __str__(self):
        return str(self._record)

    def bnd_direction_strand(self):
        """
        Parses a BND ALT string to determine its direction and strand.
         ALT  Meaning
        t[p[ piece extending to the right of p is joined after t
        t]p] reverse comp piece extending left of p is joined after t
        ]p]t piece extending to the left of p is joined before t
        [p[t reverse comp piece extending right of p is joined before t

        A direction of 'left' means the piece is anchored on the left of the breakpoint

        Note that this method assumes `self.is_bnd()`

        Returns:
            tuple: A tuple containing the direction ("left" or "right") and strand ("direct" or "complement").
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

        Returns:
            tuple: A tuple containing the chromosome (str) and position (int).
        """
        # Regular expression to match the BND format and extract chrom:pos
        match = re.search(r'[\[\]]([^\[\]:]+):(\d+)[\[\]]', self.alts[0])
        if not match:
            raise ValueError(f"Invalid BND ALT format: {self.alts[0]}")

        chrom = match.group(1)  # Extract the chromosome
        pos = int(match.group(2))  # Extract the position as an integer

        return chrom, pos

    def bnd_match(self, other, matcher=None, **_kwargs):
        """
        Build a MatchResult for bnds
        """
        if matcher is None:
            matcher = truvari.Matcher()

        def bounds(entry, pos, key):
            """
            Inflate a bnd position based on CIPOS.
            """
            start = pos - matcher.params.bnddist
            end = pos + matcher.params.bnddist

            key = 'CI' + key
            idx = 0 if key == 'POS' else 1
            if key in entry.info:  # Just CIPOS
                start -= abs(entry.info[key][idx])
                end += abs(entry.info[key][idx])

            return start, end

        ret = truvari.MatchResult()
        ret.base = self
        ret.comp = other

        # Only put start distance same chrom pos2
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

        matcher.compare_gts(ret, self, other)

        # Score is percent of allowed distance needed to find this match
        if matcher.params.bnddist > 0:
            ret.score = max(0, (1 - ((abs(ret.st_dist) + abs(ret.ed_dist)) / 2)
                                / matcher.params.bnddist) * 100)
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
            size = self.size()
            start -= size // 2
            end += size // 2
        return start, end

    def distance(self, other):
        """
        Calculate the start and end distances of the pair. Negative distances
        indicate entryA is upstream of entryB

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
        Returns None if GT is not present
        """
        if "GT" not in self.samples[sample]:
            return None
        return self.samples[sample]['GT']

    def get_ref(self):
        """
        Returns either the resolved REF or _record.ref
        """
        if self.resolved_ref:
            return self.resolved_ref
        return self._record.ref

    def get_alt(self):
        """
        Returns either the resolved ALT or _record.alt[0]
        """
        if self.resolved_alt:
            return self.resolved_alt
        return self._record.alts[0]

    def is_bnd(self):
        """
        Returns if a record is a resolved BND
        """
        return self.alleles_variant_types[1] == 'BND'

    def is_filtered(self, values=None):
        """
        Checks if entry should be filtered given the filter values provided.
        If values is None, assume that filter must have PASS or be blank '.'

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

    def is_present(self, sample=None, allow_missing=True):
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
            >>> import pysam
            >>> v = truvari.VariantFile('repo_utils/test_files/variants/input1.vcf.gz')
            >>> e = next(v)
            >>> e.is_present()
            True
        """
        gt = self.gt(sample)
        if allow_missing:
            return 1 in gt
        # Hemi...
        return truvari.get_gt(gt)[truvari.GT.HET, truvari.GT.HOM]

    def is_single_bnd(self):
        """
        Returns if a record is a single-end BND
        """
        return '.' in (self.alts[0][0], self.alts[0][-1])

    def is_resolved(self):
        """
        Checks if an SV is a sequence resolved by checking it is
        has an alt, is not symbolic, is not monref, is not BND

        :return: is_resolved
        :rtype: bool
        """
        return self.alts and \
            ('<' not in self.alts[0]) and (self.alts[0] not in (None, '*')) \
            and self.alleles_variant_types[1] != 'BND'

    def get_record(self):
        """
        Return internal pysam.VariantRecord
        """
        return self._record

    def match(self, other, matcher=None, skip_gt=False, short_circuit=False):
        """
        Build a MatchResult from comparison of two VariantRecords
        If self and other are non-bnd, calls VariantRecord.var_match,
        If self and other are bnd, calls VariantRecord.bnd_match
        Otherwise, raises a TypeError.
        If no matcher is provided, builds one from defaults
        """
        if matcher is None:
            matcher = truvari.Matcher()
        if not self.is_bnd() and not other.is_bnd():
            return self.var_match(other, matcher, skip_gt, short_circuit)
        if self.is_bnd() and other.is_bnd():
            return self.bnd_match(other, matcher)
        raise TypeError("Incompatible Variants (BND and !BND) can't be matched")

    def move_record(self, out_vcf, sample=None):
        """
        Similar to pysam.VariantRecord.translate, except this allows adding new
        FORMAT fields and optionally selects a subset of samples
        """
        ret = out_vcf.new_record(contig=self.contig, start=self.start, stop=self.stop, alleles=self.alleles,
                                 id=self.id, qual=self.qual, filter=self.filter, info=self.info)

        for k, v in self.samples[sample].items():
            ret.samples[0][k] = v
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
            >>> import pysam
            >>> v = truvari.VariantFile('repo_utils/test_files/variants/input1.vcf.gz')
            >>> a = next(v)
            >>> b = next(v)
            >>> a.recovl(b)
            0
        """
        astart, aend = self.boundaries(ins_inflate)
        bstart, bend = other.boundaries(ins_inflate)
        return truvari.reciprocal_overlap(astart, aend, bstart, bend)

    def resolve(self, ref, dup_to_ins=False):
        """
        Attempts to resolve an SV's REF/ALT sequences. Stores in self.resolved_ref, self.resolved_alt, and self.end
        """
        if ref is None or self.alts[0] in ['<CNV>', '<INS>'] or self.start > ref.get_reference_length(self.chrom):
            return False

        seq = ref.fetch(self.chrom, self.start, self.end)
        svtype = self.var_type()
        if svtype == truvari.SV.DEL:
            self.resolved_ref = seq
            self.resolved_alt = seq[0]
        elif svtype == truvari.SV.INV:
            self.resolved_ref = seq
            self.resolved_alt = seq.translate(RC)[::-1]
        elif svtype == truvari.SV.DUP and dup_to_ins:
            self.resolved_ref = seq[0]
            self.resolved_alt = seq
            self.end = self.start + 1
        else:
            return False

        return True

    def same_type(self, other, dup_to_ins=False):
        """
        Check if entryA svtype == entryB svtype

        :param `other`: Other entry to compare with
        :type `other`: :class:`truvari.VariantRecord`
        :param `dup_to_ins`: Convert DUP to INS types
        :type `dup_to_ins`: bool

        :return: True if entry SVTYPEs match
        :rtype: bool
        """
        a_type = self.var_type()
        b_type = other.var_type()
        if dup_to_ins and a_type == truvari.SV.DUP:
            a_type = truvari.SV.INS
        if dup_to_ins and b_type == truvari.SV.DUP:
            b_type = truvari.SV.INS
        return a_type == b_type

    def seqsim(self, other, roll=True):
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

        a_seq = self.ref if self.var_type() == truvari.SV.DEL else self.get_alt()
        a_seq = a_seq.upper()
        b_seq = other.ref if other.var_type() == truvari.SV.DEL else other.get_alt()
        b_seq = b_seq.upper()
        st_dist, ed_dist = self.distance(other)

        if not roll or st_dist == 0 or ed_dist == 0:
            return truvari.seqsim(a_seq, b_seq)

        if st_dist < 0:
            st_dist *= -1
        else:
            a_seq, b_seq = b_seq, a_seq

        # Return best of rolled, unrolled from both ends, and direct similarity
        # Whichever is highest is how similar these sequences can be
        return truvari.best_seqsim(a_seq, b_seq, st_dist)

    def size(self):
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
        if "SVLEN" in self.info:
            if type(self.info["SVLEN"]) in [list, tuple]:
                size = abs(self.info["SVLEN"][0])
            else:
                size = abs(self.info["SVLEN"])
        elif self.alts is not None and self.alts[0].count("<"):
            start, end = self.boundaries()
            size = end - start
        else:
            r_len = len(self.ref)
            a_len = len(self.alts[0]) if self.alts is not None else 0
            if r_len == a_len:
                if r_len == 1:
                    size = 0  # SNPs are special
                else:
                    size = r_len
            else:
                size = abs(r_len - a_len)
        return size

    def sizesim(self, other):
        """
        Calculate the size similarity and difference for two entries

        :param `other`: Other entry to compare with
        :type `other`: :class:`truvari.VariantRecord`

        :return: size similarity and size diff (A - B)
        :rtype: (float, int)

        Example
            >>> import truvari
            >>> import pysam
            >>> v = truvari.VariantFile('repo_utils/test_files/variants/input1.vcf.gz')
            >>> a = next(v)
            >>> b = next(v)
            >>> a.sizesim(b)
            (0.07142857142857142, 13)
        """
        return truvari.sizesim(self.size(), other.size())

    def var_match(self, other, matcher=None, skip_gt=False, short_circuit=False):
        """
        Build a MatchResult
        if skip_gt, don't do genotype comparison
        if short_circuit, return after first failure
        """
        ret = truvari.MatchResult()
        ret.base = self
        ret.comp = other

        ret.state = True

        if not matcher.params.typeignore and not self.same_type(other, matcher.params.dup_to_ins):
            logging.debug("%s and %s are not the same SVTYPE",
                          str(self), str(other))
            ret.state = False
            if short_circuit:
                return ret

        bstart, bend = self.boundaries()
        cstart, cend = other.boundaries()
        if not truvari.overlaps(bstart - matcher.params.refdist, bend + matcher.params.refdist, cstart, cend):
            logging.debug("%s and %s are not within REFDIST",
                          str(self), str(other))
            ret.state = False
            if short_circuit:
                return ret

        ret.sizesim, ret.sizediff = self.sizesim(other)
        if ret.sizesim < matcher.params.pctsize:
            logging.debug("%s and %s size similarity is too low (%.3f)",
                          str(self), str(other), ret.sizesim)
            ret.state = False
            if short_circuit:
                return ret

        if not skip_gt:
            matcher.compare_gts(ret, self, other)

        ret.ovlpct = self.recovl(other)
        if ret.ovlpct < matcher.params.pctovl:
            logging.debug("%s and %s overlap percent is too low (%.3f)",
                          str(self), str(other), ret.ovlpct)
            ret.state = False
            if short_circuit:
                return ret

        if matcher.params.pctseq > 0:
            ret.seqsim = self.seqsim(other, matcher.params.no_roll)
            if ret.seqsim < matcher.params.pctseq:
                logging.debug("%s and %s sequence similarity is too low (%.3ff)",
                              str(self), str(other), ret.seqsim)
                ret.state = False
                if short_circuit:
                    return ret
        else:
            ret.seqsim = 0

        ret.st_dist, ret.ed_dist = bstart - cstart, bend - cend
        ret.calc_score()

        return ret

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
        sv_alt_match = re.compile(r"\<(?P<SVTYPE>.*)\>")

        try:
            if self.alleles_variant_types[1] == 'BND':
                return truvari.SV.BND
        except (IndexError, AttributeError):
            pass

        if "SVTYPE" in self.info:
            ret_type = self.info["SVTYPE"]
            if isinstance(ret_type, (list, tuple)):
                ret_type = ret_type[0]
            return truvari.get_svtype(ret_type)

        if self.is_resolved():
            if len(self.ref) < len(self.alts[0]):
                ret_type = "INS"
            elif len(self.ref) > len(self.alts[0]):
                ret_type = "DEL"
            elif len(self.ref) == len(self.alts[0]):
                ret_type = "SNP" if len(self.ref) == 1 else "UNK"
            return truvari.get_svtype(ret_type)

        mat = sv_alt_match.match(
            self.alts[0]) if self.alts is not None else None
        if mat is not None:
            return truvari.get_svtype(mat.groupdict()["SVTYPE"].split(':')[0])

        return truvari.get_svtype("UNK")

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
        Extract entry and tree boundaries to call `coords_within`
        """
        qstart, qend = self.boundaries()
        m_ovl = tree[self.chrom].overlap(qstart, qend)
        if len(m_ovl) != 1:
            return False
        m_ovl = list(m_ovl)[0]
        end_within = self.var_type() != truvari.SV.INS

        return truvari.coords_within(qstart, qend, m_ovl.begin, m_ovl.end - 1, end_within)

    def within(self, rstart, rend):
        """
        Extract entry boundaries and type to call `coords_within`
        """
        qstart, qend = self.boundaries()
        end_within = self.var_type() != truvari.SV.INS
        return truvari.coords_within(qstart, qend, rstart, rend, end_within)
