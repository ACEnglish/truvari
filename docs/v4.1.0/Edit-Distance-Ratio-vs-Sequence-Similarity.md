By default, Truvari uses [edlib](https://github.com/Martinsos/edlib) to calculate the edit distance between two SV calls. Optionally, the [Levenshtein edit distance ratio](https://en.wikipedia.org/wiki/Levenshtein_distance) can be used to compute the `--pctsim` between two variants. These measures are different than the sequence similarity calculated by [Smith-Waterman alignment](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm).

To show this difference, consider the following two sequences.:

```
     AGATACAGGAGTACGAACAGTACAGTACGA
     |||||||||||||||*||||||||||||||
ATCACAGATACAGGAGTACGTACAGTACAGTACGA

30bp Aligned
1bp Mismatched  (96% similarity)
5bp  Left-Trimmed (~14% of the bottom sequence)
```

The code below runs swalign, Levenshtein, and edlib to compute the `--pctsim` between the two sequences.


```python
import swalign
import Levenshtein
import edlib

seq1 = "AGATACAGGAGTACGAACAGTACAGTACGA"
seq2 = "ATCACAGATACAGGAGTACGTACAGTACAGTACGA"

scoring = swalign.NucleotideScoringMatrix(2, -1)
alner = swalign.LocalAlignment(scoring, gap_penalty=-2, gap_extension_decay=0.5)
aln = alner.align(seq1, seq2)
mat_tot = aln.matches
mis_tot = aln.mismatches
denom = float(mis_tot + mat_tot)
if denom == 0:
    ident = 0
else:
    ident = mat_tot / denom
scr = edlib.align(seq1, seq2)
totlen = len(seq1) + len(seq2)

print('swalign', ident)
# swalign 0.966666666667
print('levedit', Levenshtein.ratio(seq1, seq2))
# levedit 0.892307692308
print('edlib', (totlen - scr["editDistance"]) / totlen)
# edlib 0.9076923076923077
```

Because the swalign procedure only considers the number of matches and mismatches, the `--pctsim` is higher than the edlib and Levenshtein ratio. 

If we were to account for the 5 'trimmed' bases from the Smith-Waterman alignment when calculating the `--pctsim` by counting each trimmed base as a mismatch, we would see the similarity drop to ~83%. 

[This post](https://stackoverflow.com/questions/14260126/how-python-levenshtein-ratio-is-computed) has a nice response describing exactly how the Levenshtein ratio is computed.

The Smith-Waterman alignment is much more expensive to compute compared to the Levenshtein ratio, and does not account for 'trimmed' sequence difference. 

However, edlib is the fastest comparison method and is used by default. Levenshtein can be specified with `--use-lev` in `bench` and `collapse`.