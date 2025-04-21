Introduction
------------

Truvari's comparison engine can match variants using a wide range of thresholds. However, some alleles can produce radically different variant representations. We could dramatically lower our thresholds to identify the match, but this would cause variants from unidentical alleles to be falsely matched. 

This problem is easiest to conceptualize in the case of 'split' variants: imagine a pipeline calls a single 100bp DEL that can also be represented as two 50bp DELs. To match these variants, we would need to loosen our thresholds to `--pick multi --pctsim 0.50 --pctsize 0.50`. Plus, these thresholds leave no margin for error. If the variant caller erroneously deleted an extra base to make a 101bp DEL we would have to lower our thresholds even further. These thresholds are already too low because there's plenty of distinct alleles with >= 50% homology.

How do we deal with inconsistent representations? In an ideal world, we would simply get rid of them by harmonizing the variants. This is the aim of `truvari phab` (pronounced 'fab')

`truvari phab` is designed to remove variant representation inconsistencies through harmonization. By reconstructing haplotypes from variants, running multiple-sequence alignment of the haplotypes along with the reference, and then recalling variants, we expect to remove discordance between variant representations and simplify the work required to perform variant comparison.

Alignment Methods
-----------------
By default, `phab` will make the haplotypes and use `abpoa` to perform a multiple sequence alignment between them and the reference to harmonize variants. Note that `abpoa` may be non-deterministic in some cases.  

Another available MSA algorithm is an external call `mafft` which can have its parameters customized via e.g. `--mafft-params "--maxiterate 1000"`. While `mafft` is often a more accurate alignment technique, it isn't fast. This makes an external call to mafft v7.505. Therefore, phab expects the executable to be in the environment path. Download mafft and place its executable available in the `$PATH` [mafft](https://mafft.cbrc.jp/alignment/software/). Alternatively, you can use the Truvari [Docker container](Development#docker) which already has mafft ready for use.

If you're willing to sacrifice accuracy for a huge speed increase, you can use `--align wfa`. However, this isn't an MSA technique and instead independently aligns each haplotype, which will not produce the most parsimonious set of variants.


Example
-------
As an example, we'll use Truvari's test files in `repo_utils/test_files/phab*` which were created from real data over a tandem repeat at GRCh38 chr1:26399065-26401053 and translated to a small test genome with coordinates chr1:1-1988.

* `phab_base.vcf.gz` - an 86 sample squared-off pVCF
* `phab_comp.vcf.gz` - a single sample's VCF
* `phab_ref.fa` - a subset of the GRCh38 reference

This dataset is interesting because the `HG002` sample in `phab_base.vcf.gz` uses the same sequencing experiment ([HPRC](https://github.com/human-pangenomics/HPP_Year1_Assemblies)) as the sample `syndip` in `phab_comp.vcf.gz`, but processed with a different pipeline. And as we will see, the pipeline can make all the difference.

To start, let's use `truvari bench` to see how similar the variant calls are in this region.
```bash
truvari bench --base phab_base.vcf.gz \
	--comp phab_comp.vcf.gz \
	--sizemin 0 \
	--bSample HG002 \
	--cSample syndip \
	--no-ref a \
	--output initial_bench
```
This will compare all variants greater than 1bp ( `-s 0` includes SNPs, `-s 1` would be 1bp INDELs) from the `HG002` sample to the `syndip` sample. We're also excluding any uncalled or reference homozygous sites with `--no-ref a`. The report in `initial_bench/summary.txt` shows:
```json
{
    "TP-base": 5,
    "TP-comp": 5,
    "FP": 2,
    "FN": 22,
    "precision": 0.7142857142857143,
    "recall": 0.18518518518518517,
    "f1": 0.2941176470588235,
}
```

These variants are pretty poorly matched, especially considering the `HG002` and `syndip` samples are using the same sequencing experiment. We can also inspect the `initial_bench/fn.vcf.gz` and see a lot of these discordant calls are concentrated in a 200bp window. 

Let's use `truvari phab` to harmonize the variants in this region.
```bash
truvari phab --reference phab_ref.fa \
	--region chr1:700-900 \
	--samples HG002,syndip \
	-o harmonized.vcf.gz \
        phab_base.vcf.gz phab_comp.vcf.gz
```

In our `harmonized.vcf.gz` we can see there are now only 9 variants. Let's run `truvari bench` again on the output to see how well the variants match after being harmonized.

```bash
truvari bench -b harmonized.vcf.gz \
	-c harmonized.vcf.gz \
	-s --no-ref a \
	--bSample HG002 \
	--cSample syndip \
	-o harmonized_bench/
```
Looking at `harmonized_bench/summary.txt` shows:
```json
{
    "TP-base": 8,
    "TP-comp": 8,
    "FP": 0,
    "FN": 0,
    "precision": 1.0,
    "recall": 1.0,
    "f1": 1.0,
}
```
Now there is no difference between our two sets of variants in this region.

For this variant call-set, `truvri phab` makes `truvari bench` overkill since the variants create identical haplotypes. In fact, we can benchmark simply by counting the genotypes.
```bash
$ bcftools query -f "[%GT ]\n" phab_result/output.vcf.gz | sort | uniq -c
      1 0/1 1/0
      1 1/0 0/1
      6 1/1 1/1
```
(We can ignore the phasing differences (`0/1` vs. `1/0`). These pipelines reported the parental alleles in a different order)

MSA as Merging
--------------

If you read the `truvari phab --help` , you may have noticed that multiple VCFs can be provided and there's a `--samples` parameter. This is by design so that we can harmonize the variants across any number of VCFs. Even a single, multi-sample VCF. By performing a multiple-sequence alignment across samples, we can better represent variation across a population. 

To see this in action, let's run `phab` on all 86 samples in the `repo_utils/test_files/phab_base.vcf.gz`
```bash
truvari phab -f phab_ref.fa \
	-r chr1:700-900 \
	-o msa_example.vcf.gz \ phab_base.vcf.gz
```

As a simple check, we can count the number of variants before/after `phab`:
```bash
bcftools view -r chr1:700-900 phab_base.vcf.gz | grep -vc "#"
bcftools view -r chr1:700-900 msa_example.vcf.gz | grep -vc "#"
```
The `160` original variants given to `phab` became just `60`.

Better yet, these fewer variants occur on fewer positions:
```bash
bcftools query -r chr1:700-900 -f "%POS\n" phab_base.vcf.gz | sort | uniq | wc -l
bcftools query -r chr1:700-900 -f "%POS\n" msa_example.vcf.gz | sort | uniq | wc -l
```
This returns that the variants were over `98` positions but now sit at just `16`

We can also observe changes in the allele frequency after running `phab`:
```bash
bcftools +fill-tags -r chr1:700-900 phab_base.vcf.gz | bcftools query -f "%AC\n" | sort -n | uniq -c
bcftools +fill-tags -r chr1:700-900 msa_example.vcf.gz | bcftools query -f "%AC\n" | sort -n | uniq -c
```
The allele-count (AC) shows a 15% reduction in singletons and removal of all variants with an AF > 0.50 which would have suggested the reference holds a minor allele.
```txt
  original  phab
    # AC     # AC
   39 1     33 1
   18 2      4 2
    3 3      2 3
    3 4      2 4
    2 5      1 5
         ...
    3 69     1 35
    1 89     1 40
    8 109    1 53
    1 132    1 56
    1 150    1 81
```

Limitations
-----------
* Creating and aligning haplotypes is impractical for very long sequences and maybe practically impossible for entire human chromosomes. Therefore, `truvari phab` is recommended to only be run on sub-regions.
* By giving the variants new representations, variant counts will likely change and all original variant metadata (e.g. `FORMAT/DP`) is lost. 
* `phab` should only be used on phased variants unless you really know what you're doing. While it can run on unphased variants, it isn't recommended. If regions contain unphased Hets or overlapping variants, it becomes more difficult to build an accurate haplotype from the variants, and in many cases the final constructed haplotype will not represent any sequence observed in a sample, which will create entirely false variants. So you can try out unphased variants, but proceed with caution.

Failed Jobs
-----------
Sometimes the haplotype alignment job can fail outside of python's control. For example, `abpoa` can log an error on jobs with a large number of haplotypes.
```
[SIMDMalloc] posix_memalign fail!
Size: 274877906944, Error: ENOMEM
```
Phab will monitor these processes and catch their errors. Therefore, it is possible for phab to have an exit status of 0, even though some (or even all) of the harmonization tasks did not complete. The variants inside these failed regions will not be in the output vcf. 

`--no-dedup`
------------
By default, phab will only harmonize one representation of each observed haplotype. This gives huge performance boosts and produces identical results for `--align wfa/poa`. However, `mafft` results are not identical between when deduplicating haplotypes because it considers the relative weight of sequences based on their frequency.