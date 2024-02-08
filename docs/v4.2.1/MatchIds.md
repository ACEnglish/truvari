MatchIds are used to tie base/comparison calls together in post-processing for debugging or other exploring. MatchIds have a structure of `{chunkid}.{callid}`. The chunkid is unique id per-chunk of calls. All calls sharing chunkid were within `--chunksize` distance and were compared. The callid is unique to a call in a chunk for each VCF. Because `bench` processes two VCFs (the base and comparison VCFs), the `MatchId` has two values: the first is the base variant's MatchId and the second the comparison variant's MatchId. 

For `--pick single`, the two MatchIds will be identical in the e.g. tp-base.vcf.gz and tp-comp.vcf.gz. However, for `--pick ac|multi`, it's possible to have cases such as one base variant matching to multiple comparison variants. That would give us MatchIds like:

```
# tp-base.vcf
MatchId=4.0,4.1

# tp-comp.vcf
MatchId=4.0,4.1
MatchId=4.0,4.2
```

This example tells us that the tp-comp variants are both pointing to `4.0` in tp-base. The tp-base variant has a higher match to the tp-comp `4.1` variant. 

One easy way to combine matched variants is to use `truvari vcf2df` to convert a benchmarking result to a pandas DataFrame and leverage pandas' merge operation. First, we convert the `truvari bench` result.

```bash
truvari vcf2df --info --bench-dir bench_result/ data.jl
```

Next, we combine rows of matched variants:
```python
import joblib
import pandas as pd

# Load the data
data = joblib.load("data.jl")

# Separate out the variants from the base VCF and add new columns of the base/comp ids
base = data[data['state'].isin(['tpbase', 'fn'])].copy()
base['base_id'] = base['MatchId'].apply(lambda x: x[0])
base['comp_id'] = base['MatchId'].apply(lambda x: x[1])

# Separate out the variants from the comparison VCF and add new columns of the base/comp ids
comp = data[data['state'].isin(['tp', 'fp'])].copy()
comp['base_id'] = comp['MatchId'].apply(lambda x: x[0])
comp['comp_id'] = comp['MatchId'].apply(lambda x: x[1])

# Merge the base/comparison variants
combined = pd.merge(base, comp, left_on='base_id', right_on='comp_id', suffixes=('_base', '_comp'))

# How many comp variants matched to multiple base variants?
counts1 = combined['base_id_comp'].value_counts()
print('multi-matched comp count', (counts1 != 1).sum())

# How many base variants matched to multiple comp variants?
counts2 = combined['comp_id_base'].value_counts()
print('multi-matched base count', (counts2 != 1).sum())
```

The `MatchId` is also used by `truvari collapse`. However there are two differences. First, in the main `collapse` output, the relevant INFO field is named `CollapsedId`. Second, because collapse only has a single input VCF, it is much easier to merge DataFrames. To merge collapse results kept variants with those that were removed, we again need to convert the VCFs to DataFrames:

```bash
truvari vcf2df -i kept.vcf.gz kept.jl
truvari vcf2df -i removed.vcf.gz remov.jl
```

Then we combine them:
```python
import joblib
import pandas as pd

# Load the kept variants and set the index.
kept = joblib.load("kept.jl").set_index('CollapseId')

# Load the removed variants and set the index.
remov = joblib.load("remov.jl")
remov['CollapseId'] = remov['MatchId'].apply(lambda x: x[0])
remov.set_index('CollapseId', inplace=True)

# Join the two sets of variants
result_df = kept.join(remov, how='right', rsuffix='_removed')
```