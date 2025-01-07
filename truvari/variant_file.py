"""
Truvari class wrapping pysam VariantFile
"""
import pysam

import truvari
from truvari.region_vcf_iter import region_filter


class VariantFile:
    """
    Wrapper around pysam.VariantFile with helper functions for iteration.

    .. note::
        The context manager functionality of pysam.VariantFile is not available with truvari.VariantFile.
    """

    def __init__(self, filename, *args, params=None, **kwargs):
        """
        Initialize the VariantFile wrapper.

        :param filename: Path to the VCF file to be opened.
        :type filename: str
        :param params: Parameters to apply to all VariantRecords
        :type params: `truvari.VariantParams`
        :param args: Additional positional arguments to pass to pysam.VariantFile.
        :param kwargs: Additional keyword arguments to pass to pysam.VariantFile.
        """
        self.params = params
        self._vcf = pysam.VariantFile(filename, *args, **kwargs)

    def __getattr__(self, name):
        """
        Delegate attribute access to the original VariantFile.

        :param name: Attribute name to access from the underlying pysam.VariantFile object.
        :type name: str
        :return: The requested attribute from the underlying pysam.VariantFile.
        :rtype: Any
        """
        return getattr(self._vcf, name)

    def __iter__(self):
        """
        Iterate over the `pysam.VariantFile`, wrapping entries into `truvari.VariantRecord`.

        :return: Iterator of truvari.VariantRecord objects.
        :rtype: iterator
        """
        for i in self._vcf:
            yield truvari.VariantRecord(i, self.params)

    def __next__(self):
        """
        Return the next truvari.VariantRecord in the VariantFile.

        :return: Next truvari.VariantRecord.
        :rtype: truvari.VariantRecord
        """
        return truvari.VariantRecord(next(self._vcf), self.params)

    def fetch(self, *args, **kwargs):
        """
        Fetch variants from the `pysam.VariantFile`, wrapping them into `truvari.VariantRecords`.

        :param args: Positional arguments for the pysam.VariantFile.fetch method.
        :param kwargs: Keyword arguments for the pysam.VariantFile.fetch method.
        :return: Iterator of truvari.VariantRecord objects.
        :rtype: iterator
        """

        for i in self._vcf.fetch(*args, **kwargs):
            yield truvari.VariantRecord(i, self.params)

    def fetch_bed(self, bed_fn, inside=True, with_region=False):
        """
        Fetch variants from the VCF based on regions defined in a BED file.

        :param bed_fn: Path to the BED file defining regions of interest.
        :type bed_fn: str
        :param inside: If True, fetch variants inside the regions. If False, fetch variants outside the regions.
        :type inside: bool
        :param with_region: If True, return tuples of (`truvari.VariantRecord`, region). Defaults to False.
        :type with_region: bool
        :return: Iterator of truvari.VariantRecord objects or tuples of (`truvari.VariantRecord`, region).
        :rtype: iterator
        """
        tree = truvari.read_bed_tree(bed_fn)
        return self.fetch_regions(tree, inside, with_region)


    def fetch_regions(self, tree, inside=True, with_region=False):
        """
        Fetch variants from the VCF based on regions defined in a tree of chrom:IntervalTree.

        :param tree: Tree of chrom:IntervalTree defining regions of interest.
        :type tree: dict
        :param inside: If True, fetch variants inside the regions. If False, fetch variants outside the regions.
        :type inside: bool
        :param with_region: If True, return tuples of (`truvari.VariantRecord`, region). Defaults to False.
        :type with_region: bool
        :return: Iterator of truvari.VariantRecord objects or tuples of (`truvari.VariantRecord`, region).
        :rtype: iterator
        """
        return region_filter(self, tree, inside, with_region)

    def write(self, record):
        """
        Write a `truvari.VariantRecord` to the `pysam.VariantFile`.

        :param record: The truvari.VariantRecord to be written.
        :type record: `truvari.VariantRecord`
        """
        out = record.get_record()
        if self.params and self.params.write_resolved:
            out.ref = record.get_ref()
            out.alts = (record.get_alt(),)
        self._vcf.write(out)
