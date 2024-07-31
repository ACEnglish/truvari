`stratify` is a helper utility for counting variants within bed regions which is essentially the same as running `bedtools intersect -c`. When working with benchmarking results, there are are four vcfs to count (tp-base, tp-comp, fn, fp). Instead of running bedtools four times and collating the results, `stratify` can be given a single `bench` result directory to generate the counts.

For example:
```bash
$ truvari stratify input.bed bench/
chrom    start     end   tpbase	tp	fn	fp
chr20	280300	280900	 0	0	0	0
chr20	100000	200000	 1	1	0	0
chr20	642000	642350	 1	1	2	1
```

The output from this can then be parsed to generate more details:

```python
import pandas as pd
import truvari

df = pd.read_csv("stratify.output.txt", sep='\t')

# If the input.bed didn't have a header and so we couldn't use the `--header` parameter, we need to name columns
df.columns = ['chrom', 'start', 'end', 'tpbase', 'tp', 'fn', 'fp']

# Create the precision, recall, and f1 for each row
metrics = df[["tpbase", "tp", "fn", "fp"]].apply((lambda x: truvari.performance_metrics(*x)), axis=1)

# metrics is now a DataFrame with a single column of tuples, lets separate them into columns
metrics = pd.DataFrame(metrics.to_list(), columns=["precision", "recall", "f1"])

# Extend the dataframe's columns
df = df.join(metrics)
df.head()
```
Which gives the result:
```
   chrom   start     end  tpbase  tp  fn  fp  precision  recall        f1
0  chr20  135221  239308       1   1   0   0       1.00    1.00  1.000000
1  chr20  260797  465632       3   3   3   1       0.75    0.50  0.600000
2  chr20  465866  622410       1   1   0   0       1.00    1.00  1.000000
3  chr20  623134  655257       1   1   3   1       0.50    0.25  0.333333
4  chr20  708338  732041       1   1   1   0       1.00    0.50  0.666667
```

```
usage: stratify [-h] [-o OUT] [--header] [-w] [--debug] BED VCF

Count variants per-region in vcf

positional arguments:
  BED                   Regions to process
  VCF                   Truvari bench result directory or a single VCF

optional arguments:
  -h, --help            show this help message and exit
  -o OUT, --output OUT  Output bed-like file
  --header              Input regions have header to preserve in output
  -w, --within          Only count variants contained completely within region boundaries
  --debug               Verbose logging
```