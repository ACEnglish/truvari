MatchIds are used to tie base/comparison calls together in post-processing for debugging or other exploring. MatchIds are in `bench` and `collapse` outputs. In the `collapse` main output, the field is named 'CollapseId'.

MatchIds have a structure of `{chunkid}.{baseid}.{compid}`.

Where:
 - chunkid - Unique id per-chunk of calls. All calls sharing chunkid were within `--chunksize` distance and were
   compared.
 - baseid - The number of the base call from this chunk
 - compid - The number of the comp call from this chunk

Doing basic pairing with MatchIds when results didn't use `--pick multi` is straight forward. For example, to tie tp-base.vcf and tp-comp.vcf calls together, they will have identical MatchIds.

MatchIds with `_` for the baseid or compid are calls that were in a chunk without a counterpart (e.g. `5._.0` is the
first comparison call from the fifth chunk, which had no base calls.)

When using `--pick multi` there can be instances such as:
```
VCF         MatchId
tp-base     3.0.0
tp-comp     3.0.0
tp-comp     3.0.1
```
Here, the 0th base/comp call match to one another. Additionally, the 1st comp call `3.0.1` also matches to the 0th base call `tp-base 3.0.0`.

This makes tieing false positives/negatives to their closest matching counter part a little tricker.
For example, if we have a false negative with MatchId `3.1.0`, we'll need to look in both the `fp.vcf` and `tp-comp.vcf` for
`3.*.0`. This search will find the 0th comp call in the third chunk which may have matched to a different baseid than 1.