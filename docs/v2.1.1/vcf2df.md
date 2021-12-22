# vcf2df

After you've created your benchmarking results with `truvari bench`, you'll often want to plot different views of your results. We enjoy using [pandas](https://pandas.pydata.org/)/[seaborn](https://seaborn.pydata.org/) for python plotting, so we've made the command `truvari vcf2df`. This will take an output Truvari directory and convert each vcf entry into a row of a pandas DataFrame. The resulting DataFrame will always have the columns:
* id : VCF column ID
* svtype : SVTYPE as determined by `truvari.entry_variant_type`
* svlen : SVLEN as determined by `truvari.entry_size`
* szbin : SVLEN's size bin as determined by `truvari.get_sizebin`
* qual : VCF column QUAL
* filter : VCF column FILTER
* is_pass : boolean of if the filter is empty or PASS
* state : The truvari state assigned to the variant 
    * tpbase : Parsed from the tp-base.vcf
    * tp : Parsed from the tp-call.vcf
    * fp : Parsed from the fp.vcf
    * fn : Parsed from the fn.vcf

Additionally, the DataFrame is assigned an index of a description of the variant `CHROM:START-STOP.ALT`. This enables operations such as `pandas.DataFrame.join` to assist creating a union between sets of variants. 

Optionally, `vcf2df` can attempt to pull `INFO` and `FORMAT` fields from the VCF and put each field into the DataFrame as a new column. For FORMAT fields, the VCF header definition's `Number` is considered and multiple columns may be added. For example, the `AD` field, typically holding Allele Depth has `Number=A`, indicating that there will be one value for each allele. Truvari assumes that all VCFs are single-allele per line, so there are only 2 alleles described per-VCF line, the reference and alternate allele. Therefore, two columns are added to the DataFrame, `AD_ref` and `AD_alt` corresponding to the 0th and 1st values from the AD field's list of values. Similarity, for PL (genotype likelihood) with `Number=G`, there's three values and columns are created named `PL_ref`, `PL_het`, `PL_hom`. 

The created DataFrame is saved into a joblib file, which can then be plotted as simply as:
```python
import joblib
import seaborn as sb
import matplotlib.pyplot as plt

data = joblib.load("test.jl")
p = sb.countplot(data=data[data["state"] == 'tp'], x="szbin", hue="svtype", hue_order=["DEL", "INS"])
plt.xticks(rotation=45, ha='right')
p.set(title="True Positives by svtype and szbin")
```
![](https://github.com/spiralgenetics/truvari/blob/develop/imgs/truv2df_example.png)

This enables concatenation of Truvari results across multiple benchmarking experiments for advanced comparison. For example, imagine there's multiple parameters used for SV discovery over multiple samples. After running `truvari bench` on each of the results with the output directories named to `params/sample/` and each converted to DataFrames with `truvari vcf2df`, we can expand/concatenate the saved joblib DataFrames with:

```
import glob
import joblib
import pandas as pd

files = glob.glob("*/*/data.jl")
dfs = []
for f in files:
    params, sample, frame = f.split('/')
    d = joblib.load(f)
    d["params"] = params
    d["sample"] = sample
    dfs.append(d)
df = pd.concat(dfs)
joblib.dump(df, "results.jl")
```

Since these DataFrames may be useful outside of Truvari results, by default a single VCF file is converted to a DataFrame where `--bench-dir` needs to be specified for a directory.

```bash
usage: vcf2df [-h] [-b] [-i] [-f] [-S] [--debug] VCF JL

Takes a vcf and creates a data frame. Can parse a bench output directory

positional arguments:
  VCF                   VCF to parse
  JL                    Output joblib to save

optional arguments:
  -h, --help            show this help message and exit
  -b, --bench-dir       Input is a truvari bench directory
  -i, --info            Attempt to put the INFO fields into the dataframe
  -f, --format          Attempt to put the FORMAT fileds into the dataframe
  -S, --skip-compression
                        Skip the attempt to optimize the dataframe's size
  --debug               Verbose logging
```