A frequent application of comparing SVs is to perform a 'bakeoff' of performance
between two SV programs against a single set of base calls.

Beyond looking at the Truvari results/report, you may like to investigate what calls
are different between the programs.

Below is a set of scripts that may help you generate those results. For our examples,
we'll be comparing arbitrary programs Asvs and Bsvs aginst base calls Gsvs.

*_Note_* - This assumes that each record in Gsvs has a unique ID in the vcf.

Generate the Truvari report for Asvs and Bsvs
=============================================

```bash
truvari bench -b Gsvs.vcf.gz -c Asvs.vcf.gz -o cmp_A/ ...
truvari bench -b Gsvs.vcf.gz -c Bsvs.vcf.gz -o cmp_B/ ...
```
Consistency
===========
The simplest way to compare the programs is to get the intersection of TPbase calls from the two reports.
```bash
truvari consistency cmp_A/tp-base.vcf cmp_B/tp-base.vcf
```
See [[consistency wiki|consistency]] for details on the report created.

Below are older notes to manually create a similar report to what one can make using `truvari consistency`

Combine the TPs within each report
==================================

```bash
cd cmp_A/
paste <(grep -v "#" tp-base.vcf) <(grep -v "#" tp-comp.vcf) > combined_tps.txt
cd ../cmp_B/
paste <(grep -v "#" tp-base.vcf) <(grep -v "#" tp-comp.vcf) > combined_tps.txt
```

Grab the FNs missed by only one program
=======================================

```bash
(grep -v "#" cmp_A/fn.vcf && grep -v "#" cmp_B/fn.vcf) | cut -f3 | sort | uniq -c | grep "^ *1 " | cut -f2- -d1 > missed_names.txt
```

Pull the TP sets' difference
============================

```bash
cat missed_names.txt | xargs -I {} grep -w {} cmp_A/combined_tps.txt > missed_by_B.txt
cat missed_names.txt | xargs -I {} grep -w {} cmp_B/combined_tps.txt > missed_by_A.txt
```

To look at the base-calls that Bsvs found, but Asvs didn't, run `cut -f1-12 missed_by_A.txt`.

To look at the Asvs that Bsvs didn't find, run `cut -f13- missed_by_B.txt`.

Shared FPs between the programs
===============================

All of the work above has been about how to analyze the TruePositives. If you'd like to see which calls are shared between Asvs and Bsvs that aren't in Gsvs, simply run Truvari again.

```bash
bgzip cmp_A/fp.vcf && tabix -p vcf cmp_A/fp.vcf.gz
bgzip cmp_B/fp.vcf && tabix -p vcf cmp_B/fp.vcf.gz
truvari bench -b cmp_A/fp.vcf.gz -c cmp_B/fp.vcf.gz -o shared_fps ...
```