import pysam
import truvari
import hashlib

RC = str.maketrans("ATCG", "TAGC")

class VariantRecord:
    def __init__(self, record):
        self._record = record

    # Delegate attribute access to the original record
    def __getattr__(self, name):
        return getattr(self._record, name)
    
    def gt(self, sample=0):
        if "GT" not in self.samples[sample]:
            return None
        return self.samples[sample]['GT']

    def is_bnd(self):
        return self.alleles_variant_types[1] == 'BND'
    
    def is_single_bnd(self):
        return '.' in (self.alts[0][0], self.alts[0][-1])

    def is_monrefstar(self):
        return not self.alts or self.alts[0] in (None, '*')

    def is_multi(self):
        return len(self.alts) > 1

    def svtype(self):
        return truvari.entry_variant_type(self)

    def size(self):
        return truvari.entry_size(self)

    def boundaries(self, ins_inflate=False):
        return truvari.entry_boundaries(self, ins_inflate)
    
    def distance(self, other):
        return truvari.entry_distance(self, other)

    def gt_comp(self, other, sampleA=None, sampleB=None):
        return truvari.entry_gt_comp(self, other)

    def is_filtred(self, values=None):
        return truvari.entry_is_filtered(self, values)

    def is_present(self, sample=None, allow_missing=True):
        return truvari.entry_is_present(self, sample, allow_missing)

    def seqsim(self, other, roll=True):
        return truvari.entry_seq_similarity(self, other, roll)

    def recovl(self, other, ins_inflate=True):
        return truvari.entry_reciprocal_overlap(self, other, ins_inflate)

    def same_type(self, other, dup_to_ins=False):
        return truvari.entry_same_variant_type(self, other, dup_to_ins)

    def sizesim(self, other):
        return truvari.entry_size_similarity(self, other)

    def to_hash(self, hasher=hashlib.sha1):
        return truvari.entry_to_hash(self, hasher)

    def to_key(self, prefix="", bounds=False):
        return truvari.entry_to_key(self, prefix, bounds)

    def within_tree(self, tree):
        return truvari.entry_within_tree(self, tree)

    def within(self, rstart, rend):
        return truvari.entry_within(self, rstart, rend)

    def is_resolved(self):
        return truvari.entry_resolved(self)
    
    def __str__(self):
        return str(self._record)
    
    def resolve(self, ref, dup_to_ins=False):
        """
        Attempts to resolve an SV's REF/ALT sequences
        """
        entry = self._record
        if ref is None or entry.alts[0] in ['<CNV>', '<INS>'] or entry.start > ref.get_reference_length(entry.chrom):
            return False

        # BNDs which describe a deletion can be resolved
        if entry.alleles_variant_types[1] == 'BND':
            if "SVTYPE" in entry.info and entry.info['SVTYPE'] == 'DEL':
                entry.alts = ['<DEL>']
            else:
                return False

        seq = ref.fetch(entry.chrom, entry.start, entry.stop)
        svtype = self.svtype()
        if svtype == truvari.SV.DEL:
            entry.ref = seq
            entry.alts = [seq[0]]
        elif svtype == truvari.SV.INV:
            entry.ref = seq
            entry.alts = [seq.translate(RC)[::-1]]
        elif svtype == truvari.SV.DUP and dup_to_ins:
            entry.ref = seq[0]
            entry.alts = [seq]
            entry.stop = entry.start + 1
        else:
            return False

        return True


if __name__ == '__main__':
    v = pysam.VariantFile("../repo_utils/test_files/variants/input1.vcf.gz")
    e = VariantRecord(next(v))
    print(str(e))
    print(e.ref)
    e.ref = 'something'
    print(str(e))
    print(e.ref)
