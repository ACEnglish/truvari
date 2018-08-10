import pyfaidx

ref = pyfaidx.Fasta("reference_lambda.fasta")
print dir(ref)

a1_chrom, a1_start, a1_end, a1_seq = "lambda", 500, 1000, ""
a2_chrom, a2_start, a2_end, a2_seq = "lambda", 400, 900, ""

start = min(a1_start, a2_start)
end = max(a1_end, a2_end)

print ref[a1_chrom][start:a1_start].seq
hap1_seq = ref.get_seq(a1_chrom, start + 1, a1_start).seq
hap1_seq += a1_seq
hap1_seq += ref.get_seq(a1_chrom, a1_end + 1, end).seq

hap2_seq = ref.get_seq(a2_chrom, start + 1, a2_start).seq
hap2_seq += a2_seq
hap2_seq += ref.get_seq(a2_chrom, a2_end + 1, end).seq


import Levenshtein

print Levenshtein.ratio(hap1_seq, hap2_seq)

fh = open("reference_lambda.fasta")
fh.readline()
seq = fh.readline().strip()

# no plus one
h1_oseq = seq[start:a1_start] + a1_seq + seq[a1_end:end]
h2_oseq = seq[start:a2_start] + a2_seq + seq[a2_end:end]
print Levenshtein.ratio(h1_oseq, h2_oseq)

print hap1_seq
print h1_oseq
print hap2_seq
print h2_oseq
