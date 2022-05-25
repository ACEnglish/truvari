As part of the GIAB analysis, tandem repeat (TR) stratifications for false positives (FP), false negatives (FN), and true positives (TP) were analyzed. See [this discussion](https://groups.google.com/d/msg/giab-analysis-team/tAtVBm9Fdrw/I2qeazh8AwAJ) for details.

The stratifications.sh script automates the procedure above by taking a Truvari bench output directory as the first argument and TR repeat annotations (found in the link above) to create the following files in the output directory.

Running this script requires `bedtools`, `vcf-sort`, and `bgzip` to be in your environment.


| FileName            | Description                           |
|---------------------|---------------------------------------|
| fn_ins_nonTR.vcf.gz | FN insertions not found in TR regions |
| fn_ins_TR.vcf.gz    | FN insertions found in TR regions     |
| fp_ins_nonTR.vcf.gz | FP insertions not found in TR regions |
| fp_ins_TR.vcf.gz    | FP insertions found in TR regions     |
| fn_del_nonTR.vcf.gz | FN insertions not found in TR regions |
| fn_del_TR.vcf.gz    | FN insertions found in TR regions     |
| fp_del_nonTR.vcf.gz | FP insertions not found in TR regions |
| fp_del_TR.vcf.gz    | FP insertions found in TR regions     |