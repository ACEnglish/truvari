import pysam
import truvari
import time

a = time.time()
v = pysam.VariantFile("sample.vcf.gz")
x = truvari.Matcher()

up_record = next(v)
for entry in v:
    if x.filter_call(entry):
        continue
    print("size=%s type=%s bound=%s dist=%s ovl=%s szsim=%s gtcomp=%s pres=%s filt=%s sim=%s" % (truvari.entry_size(entry),
                                        truvari.entry_variant_type(entry),
                                        truvari.entry_boundaries(entry, False),
                                        truvari.entry_distance(up_record, entry),
                                        round(truvari.entry_reciprocal_overlap(up_record, entry), 7),
                                        tuple(round(_, 7) for _ in truvari.entry_size_similarity(up_record, entry)),
                                        truvari.entry_gt_comp(up_record, entry, 0, 0),
                                        truvari.entry_is_present(entry),
                                        truvari.entry_is_filtered(entry),
                                        round(truvari.entry_seq_similarity(up_record, entry), 7),
                                        ))
    x.build_match(up_record, entry)
    up_record = entry

b = time.time()
print(b - a)

