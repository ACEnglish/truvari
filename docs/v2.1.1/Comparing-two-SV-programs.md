A frequent application of comparing SVs is to perform a 'bakeoff' of performance
between two SV programs against a single set of base calls.

Beyond looking at the Truvari results/report, you may like to investigate what calls
are different between the programs.

Below is a set of scripts that may help you generate those results. For our examples,
we'll be comparing arbitrary programs Asvs and Bsvs aginst base calls Gsvs.

*_Note_* - This assumes that each record in Gsvs has a unique ID in the vcf.

Generate the Truvari report for Asvs and Bsvs.
----------------------------------------------

```bash
./truvari.py -b Gsvs.vcf.gz -c Asvs.vcf.gz -o cmp_A/
./truvari.py -b Gsvs.vcf.gz -c Bsvs.vcf.gz -o cmp_B/
```

Combine the TPs within each report.
-----------------------------------

```bash
cd cmp_A/
paste <(grep -v "#" tp-base.vcf) <(grep -v "#" tp-call.vcf) > combined_tps.txt
cd ../cmp_B/
paste <(grep -v "#" tp-base.vcf) <(grep -v "#" tp-call.vcf) > combined_tps.txt
```

Grab the FNs missed by only one program.
----------------------------------------

```bash
(grep -v "#" cmp_A/fn.vcf && grep -v "#" cmp_B/fn.vcf) | cut -f3 | sort | uniq -c | grep "^ *1 " | cut -f2- -d1 > missed_names.txt
```

Pull the TP sets' difference.
-----------------------------

```bash
cat missed_names.txt | xargs -I {} grep -w {} cmp_A/combined_tps.txt > missed_by_B.txt
cat missed_names.txt | xargs -I {} grep -w {} cmp_B/combined_tps.txt > missed_by_A.txt
```

To look at the base-calls that Bsvs found, but Asvs didn't, run `cut -f1-12 missed_by_A.txt`.

To look at the Asvs that Bsvs didn't find, run `cut -f13- missed_by_B.txt`.

Calculate the overlap.
----------------------

One may wish for summary numbers of how many calls are shared/unique between the two programs.
Truvari has a program to help. See https://github.com/spiralgenetics/truvari/wiki/Consistency-report for details.

Shared FPs between the programs
-------------------------------

All of the work above has been about how to analyze the TruePositives. If you'd like to see which calls are shared between Asvs and Bsvs that aren't in Gsvs, simply run Truvari again.

```bash
bgzip cmp_B/fp.vcf && tabix -p vcf cmp_B/fp.vcf.gz
./truvari.py -b cmp_A/fp.vcf -c cmp_B/fp.vcf.gz -o shared_fps
```