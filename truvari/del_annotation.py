"""
Holds all the information needed to annotate a VCF from a single source
"""
import gzip
from collections import defaultdict

import pysam
import truvari.comparisons as comp
from intervaltree import IntervalTree

class AnnoData:
    """
    Holds data that can be put into an annotation
    I reckon I need to figure out a way to read the header information for these things...
    """
    return None

class Entry:
    """
    Generic entry object for use by annotators
    """
    def __init__(self, chrom, start, stop, event_type=None, event_size=None, anno_data=None):
        """
        We annotate over ranges, then you can use kwargs for each of
        """
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.anno_data = anno_data

class BEDReader:
    """
    Read bed file for use by the annotation engine
    """
    def __init__(self, fn, vcf_header=None):
        self.filename = fn
        # if vcf_header is none, then we can do stuff to try and predict what it would be?
        # Or should I require that they are a thing provided to create-anno


    def __iter__(self):
        """ iterate all the entries """
        with gzip.GzipFile(self.filename) as fh:
            header = fh.readline().decode()
            if not header.startswith("#"):
                raise Exception("Expected header line in BED file %s" % self.filename)
            header = header.strip().split('\t')

            for line in fh:
                line = line.decode()
                data = line.strip().split('\t')
                chrom, start, stop = data[:3]
                start = int(start)
                stop = int(stop)
                ann = {}
                for key, val in zip(header[3:], line[3:]):
                    ann[key] = val
                yield Entry(chrom, start, stop, ann)

"""
What are all the properties I will need 'standardized' into an Entry that will be used for comparison

I can pretty much look at everything that Truvari does with a VCF entry and replicate that
Then I can fill in sub-pieces of that with less informed pieces.
The hard part is how do I extract the information?

I could make an Annotator per-source and assume I have a vcf entry to annotate.
So the Annotator will just need to have the generic ".annotate" and then it does the work.
And I can re-use code for things like BED annotations which already exist.
And for specialized things like 1kg, I can simply tie the code to it with config and headers
And those can also be reused.
So instead of data-types, I can work on annotation sources
"""

class Annotation():
    def __init__(self, anno_file, match_ops):
        self.anno_file = anno_file
        self.match_ops = match_ops
        self.trees = defaultdict(IntervalTree)
        for entry in self.anno_file:
            self.trees[entry.chrom].addi(entry.start, entry.stop, entry)

    def annotate(self, entry):
        """
        Given an pyvcf Variant entry do the matching
        """
        #for anno_entry in self.trees[entry.chrom][entry.start:entry.stop]:
        for anno_entry in self.trees[entry.chrom].overlap(entry.start, entry.stop):
            add_annot = True
            for match in self.match_ops:
                if not match(entry, anno_entry):
                    add_annot = False
                    break
            if add_annot:
                print("would annotate %s with %s" % (entry, anno_entry.data.anno_data))
                #entry.add_annotation(anno_entry)

    def save(self, fn):
        """
        Save this as a pickle so we can just quickly use it later
        """
        pass

def test():
    """
    Quick file test
    """
    bed = BEDReader("/home/english/truvari/tests/anno/rep.bed.gz")
    match_ops = [reciprocal_overlap(.70)]

    anno = Annotation(bed, match_ops)
    for entry in pysam.VariantFile("/home/english/truvari/tests/anno/vars.vcf.gz"):
        if "SVLEN" in entry.info:
            anno.annotate(entry)

if __name__ == '__main__':
    test()
