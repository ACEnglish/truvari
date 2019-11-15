"""
Holds all the information needed to annotate a VCF from a single source
"""
import gzip
from collections import defaultdict

import truvari.comparisons as comp

from intervaltree import IntervalTree
# This holds a generic IntervalTree that has intersection methods so that we can get back the metadata

# Then we'll pickle these things up and load them on the fly...
# Operations 
#   add_intervals(chrom, start, end, **kwargs)
#      where kwargs are key of the INFO name and val is header information for the vcf,, so I should make those
#   add_match_thresholds([list, of methods, that, an entry would need to pass])
# so we need to clean up all the comparisons so that they can operate on a vcf entry, and an interval?
#   annotate(entry):
#       I'll look through my interval tree with the match thresholds and then say we got it...
#       

#GFF, bed, or vcf readers
# Need to wrap/standardize them so I can consume all of them by annotation

class AnnoData:
    """
    Holds data that can be put into an annotation
    I reckon I need to figure out a way to read the header information for these things...
    """
    pass

class Entry:
    """
    Generic entry object for use by annotators
    """
    def __init__(self, chrom, start, end, anno_data):
        """
        We annotate over ranges, then you can use kwargs for each of 
        """
        self.chrom = chrom
        self.start = start
        self.end = end
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
                chrom, start, end = data[:3]
                start = int(start)
                end = int(end)
                ann = {}
                for key, val in zip(header[3:], line[3:]):
                    ann[key] = val
                yield Entry(chrom, start, end, ann)

class Annotation():
    def __init__(self, anno_file, match_ops):
        self.anno_file = anno_file
        self.match_ops = match_ops
        self.trees = defaultdict(IntervalTree)
        for entry in self.anno_file:
            self.trees[entry.chrom].addi(entry.start, entry.end, entry)

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

# Match operations - wrappers, function wrapper whatevers for the core operations
# And you can make custom ones, too, I reckno
def reciprocal_overlap(min_pct):
    """
    creates a reciprocal overlap rule for matching two entries. Returns a method that can be used as a match operator
    """
    def ret_fun(entry1, entry2):
        return comp.get_rec_ovl(entry1.start, entry1.stop, entry2.begin, entry2.end) >= min_pct
    return ret_fun

def test():
    """
    Quick file test
    """
    bed = BEDReader("/home/english/truvari/tests/anno/rep.bed.gz")
    match_ops = [reciprocal_overlap(.70)]
    
    anno = Annotation(bed, match_ops)
    import pysam
    for entry in pysam.VariantFile("/home/english/truvari/tests/anno/vars.vcf.gz"):
        if "SVLEN" in entry.info:
            anno.annotate(entry)

if __name__ == '__main__':
    test()
