import re
import pysam
import truvari
import hashlib

RC = str.maketrans("ATCG", "TAGC")


class VariantRecord:
    """
    Wrapper around pysam.VariantRecords with helper functions of variant properties and basic comparisons
    """

    def __init__(self, record):
        """
        Initialize with just the internal record
        """
        self._record = record

    def __getattr__(self, name):
        """
        Delegate attribute access to the original record
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
            raise ValueError(f"Invalid BND ALT format: {bnd}")

        chrom = match.group(1)  # Extract the chromosome
        pos = int(match.group(2))  # Extract the position as an integer

        return chrom, pos

    def boundaries(self, ins_inflate=False):
        """
        Get the start and end of an entry

        :param `entry`: entry to get bounds
        :type `entry`: :class:`pysam.VariantRecord`
        :param `ins_inflate`: inflate INS boundaries by sv length
        :type `ins_inflate`: bool, optional

        :return: start/end
        :rtype: tuple (int, int)
        """
        start = self.start
        end = self.stop
        if ins_inflate and self.svtype() == truvari.SV.INS:
            size = self.size()
            start -= size // 2
            end += size // 2
        return start, end

    def distance(self, other):
        """
        Calculate the start and end distances of the pair. Negative distances
        indicate entryA is upstream of entryB

        :param `entryA`: first entry
        :type `entryA`: :class:`pysam.VariantRecord`
        :param `entryB`: second entry
        :type `entryB`: :class:`pysam.VariantRecord`

        :return: starts and ends distance
        :rtype: tuple (int, int)
        """
        astart, aend = entryA.boundaries()
        bstart, bend = entryB.boundaries()
        return astart - bstart, aend - bend

    def gt(self, sample=0):
        """
        Returns the raw GT field from a record for a specific sample.
        Returns None if GT is not present
        """
        if "GT" not in self.samples[sample]:
            return None
        return self.samples[sample]['GT']

    def is_bnd(self):
        """
        Returns if a record is a resolved BND
        """
        return self.alleles_variant_types[1] == 'BND'

    def is_filtered(self, values=None):
        """
        Checks if entry should be filtered given the filter values provided.
        If values is None, assume that filter must have PASS or be blank '.'

        :param `entry`: entry to check
        :type `entry`: :class:`pysam.VariantRecord`
        :param `values`: set of filter values for intersection
        :type `values`: set, optional

        :return: True if entry should be filtered
        :rtype: bool

        Example
            >>> import truvari
            >>> import pysam
            >>> v = pysam.VariantFile('repo_utils/test_files/variants/input1.vcf.gz')
            >>> e = pysam.VariantRecord(next(v)
            >>> e.is_filtered() # PASS shouldn't be filtered
            False
            >>> e = pysam.VariantRecord(next(v)
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
        return truvari.entry_is_present(self, sample, allow_missing)

    def is_single_bnd(self):
        """
        Returns if a record is a single-end BND
        """
        return '.' in (self.alts[0][0], self.alts[0][-1])

    def is_resolved(self):
        """
        Checks if an SV is a sequence resolved by checking it is
        has an alt, is not symbolic, is not monref, is not BND

        :param `entry`: entry to check
        :type `entry`: :class:`pysam.VariantRecord`

        :return: is_resolved
        :rtype: bool
        """
        return self.alts and \
            ('<' not in self.alts[0]) and not (self.alts[0] in (None, '*')) \
            and self.alleles_variant_types[1] != 'BND'

    def recovl(self, other, ins_inflate=True):
        return truvari.entry_reciprocal_overlap(self, other, ins_inflate)

    def resolve(self, ref, dup_to_ins=False):
        """
        Attempts to resolve an SV's REF/ALT sequences
        Note: This will 'replace' the truvari.VariantRecord.ref and truvari.VariantRecord.alt. To access the original
        pysam.VariantRecord attributes, use the `truvari.._record.alt`. Also, the original attributes are written to
        the any output
        """
        if ref is None or self.alts[0] in ['<CNV>', '<INS>'] or self.start > ref.get_reference_length(self.chrom):
            return False

        # BNDs which describe a deletion can be resolved
        if self.alleles_variant_types[1] == 'BND':
            if "SVTYPE" in self.info and self.info['SVTYPE'] == 'DEL':
                self.alts = ['<DEL>']
            else:
                return False

        seq = ref.fetch(self.chrom, self.start, self.stop)
        svtype = self.svtype()
        if svtype == truvari.SV.DEL:
            self.ref = seq
            self.alts = [seq[0]]
        elif svtype == truvari.SV.INV:
            self.ref = seq
            self.alts = [seq.translate(RC)[::-1]]
        elif svtype == truvari.SV.DUP and dup_to_ins:
            self.ref = seq[0]
            self.alts = [seq]
            self.stop = self.start + 1
        else:
            return False

        return True

    def same_type(self, other, dup_to_ins=False):
        return truvari.entry_same_variant_type(self, other, dup_to_ins)

    def seqsim(self, other, roll=True):
        return truvari.entry_seq_similarity(self, other, roll)

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

        :param `entry`: entry to look at
        :type `entry`: :class:`pysam.VariantRecord`

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
        return truvari.entry_size_similarity(self, other)

    def svtype(self):
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

        :param `entry`:
        :type `entry`: :class:`pysam.VariantRecord`

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
        return truvari.entry_to_hash(self, hasher)

    def to_key(self, prefix="", bounds=False):
        return truvari.entry_to_key(self, prefix, bounds)

    def within_tree(self, tree):
        return truvari.entry_within_tree(self, tree)

    def within(self, rstart, rend):
        return truvari.entry_within(self, rstart, rend)
