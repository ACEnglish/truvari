
Truvari annotations:
* [gcpct](anno#truvari-anno-gcpct) - GC Percent
* [gtcnt](anno#truvari-anno-gtcnt) - Genotype Counts
* [trf](anno#truvari-anno-trf) - Tandem Repeats
* [grm](anno#truvari-anno-grm) - Mappability
* [repmask](anno#truvari-anno-repmask) - Repeats
* [remap](anno#truvari-anno-remap) - Allele Remapping
* [hompct](anno#truvari-anno-hompct) - Homozygous Percent
* [numneigh](anno#truvari-anno-numneigh) - Number of Neighbors
* [svinfo](anno#truvari-anno-svinfo) - SVINFO Fields
* [bpovl](anno#truvari-anno-bpovl) - Annotation Intersection
* [density](anno#truvari-anno-density) - Call Density
* [dpcnt](anno#truvari-anno-dpcnt) - Depth (DP) and Alt-Depth (AD) Counts
* [lcr](anno#truvari-anno-lcr) - Low-complexity Regions
* [grpaf](anno#truvari-anno-grpaf) - Sample Group Allele-Frequency Annotations

# truvari anno gcpct

This will add an INFO tag `GCPCT` to each element in a VCF of the GC percent of the call's sequence.

For deletions, this is the GC percent of the reference range of the call. For insertions, the ALT sequence is analyzed.
```
usage: gcpct [-h] [-o OUTPUT] -r REFERENCE [input]

Annotates GC content of SVs

positional arguments:
  input                 VCF to annotate (stdin)

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output filename (stdout)
  -r REFERENCE, --reference REFERENCE
                        Reference fasta
```

# truvari anno gtcnt
This will add an INFO tag `GTCNT` to each element in a VCF with the count of genotypes found across all samples. The value is a list of Counts of genotypes for the allele across all samples (UNK, REF, HET, HOM). This is most useful for pVCFs.

```
usage: gtcnt [-h] [-o OUTPUT] [input]

Annotates GTCounts of alleles

positional arguments:
  input                 VCF to annotate (stdin)

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output filename (stdout)
```

# truvari anno trf
Adds a tandem repeat annotation to sequence resolved Insertion/Deletion variants a VCF.

### Annotations added
| Field Name | Description                                                 |
|------------|-------------------------------------------------------------|
| TRF        | Entry hits a tandem repeat region                           |
| TRFdiff    | ALT TR copy difference from reference                       |
| TRFrepeat  | Repeat motif                                                |
| TRFovl     | Percent of ALT covered by TRF annotation                    |
| TRFstart   | tart position of discovered repeat                          |
| TRFend     | End position of discovered repeat                           |
| TRFperiod  | eriod size of the repeat                                    |
| TRFcopies  | Number of copies aligned with the consensus pattern         |
| TRFscore   | Alignment score                                             |
| TRFentropy | Entropy measure                                             |
| TRFsim     | Similarity of ALT sequence to generated motif faux sequence |

### Details
Given a set of tandem repeat regions and a VCF, this annotate the tandem repeat motif and copy number change of insertions and deletions as expansions/contraction. The expected input catalog of tandem repeats is from a subset of columns in the Adotto TR catalog ([link](https://github.com/ACEnglish/adotto/blob/main/regions/DataDescription.md)). This file can be formatted for `truvari anno trf` via:
```bash
zcat adotto_TRregions_v1.1.bed.gz | cut -f1-3,18 | bgzip > anno.trf.bed.gz
tabix anno.trf.bed.gz
```
For deletions, the tool simply finds the motif annotation with the highest overlap over the variant's boundaries. It then removes that sequence from the reference and calculates how many copies of the motif are removed with the formula `round(-(ovl_pct * svlen) / anno["period"], 1)`. If a deletion overlaps multiple motifs, the highest scoring motif is chosen based on higher reciprocal overlap percent first and TRF score second (see [code](https://github.com/ACEnglish/truvari/blob/2219f52850252c18dcd8c679da6644bb1cee5b68/truvari/annotations/trf.py#L29)]. 

For insertions, by default the tool first tries to estimate which motif is contained in the alternate sequence. For each overlapping annotation, the copy number difference of the motif is calculated via `copy_diff = len(entry.alts[0][1:]) / anno["period"]`. Next, a 'feaux sequence' is made from `copy_diff` number of the motif. If the sequence is above the `--motif-similarity` with the insertion sequence, that is considered the insertion's motif. If no estimate is above the `--motif-similarity`, the insertion is incorporated into the reference and TRF is run. If the discovered TRF hits match a motif in the tandem repeat regions file, that annotation is used. If the highest scoring TRF hit doesn't match the tandem repeats region file, the nominally de novo annotation is added to the insertion's vcf entry.

```
usage: trf [-h] -i INPUT [-o OUTPUT] [-e EXECUTABLE] [-T TRF_PARAMS] -r REPEATS -f REFERENCE [-s MOTIF_SIMILARITY]
           [-m MIN_LENGTH] [-R] [--no-estimate] [-C CHUNK_SIZE] [-t THREADS] [--debug]

Intersect vcf with reference tandem repeats and annotate
variants with the best fitting repeat motif and its copy number
relative to the reference

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        VCF to annotate
  -o OUTPUT, --output OUTPUT
                        Output filename (stdout)
  -e EXECUTABLE, --executable EXECUTABLE
                        Path to tandem repeat finder (trf409.linux64)
  -T TRF_PARAMS, --trf-params TRF_PARAMS
                        Default parameters to send to trf (3 7 7 80 5 40 500 -h -ngs)
  -r REPEATS, --repeats REPEATS
                        Reference repeat annotations
  -f REFERENCE, --reference REFERENCE
                        Reference fasta file
  -s MOTIF_SIMILARITY, --motif-similarity MOTIF_SIMILARITY
                        Motif similarity threshold (0.9)
  -m MIN_LENGTH, --min-length MIN_LENGTH
                        Minimum size of entry to annotate (50)
  -R, --regions-only    Only write variants within --repeats regions (False)
  --no-estimate         Skip INS estimation procedure and run everything through TRF. (False)
  -C CHUNK_SIZE, --chunk-size CHUNK_SIZE
                        Size (in mbs) of reference chunks for parallelization (5)
  -t THREADS, --threads THREADS
                        Number of threads to use (1)
  --debug               Verbose logging
```

# truvari anno grm

For every SV, we create a kmer over the the upstream and downstream reference and alternate breakpoints. 
We then remap that kmer to the reference genome and report alignment information.
This does not alter the VCF traditional annotations, but instead will create a pandas 
DataFrame and save it to a joblib object.

There are four queries made per-SV. For both reference (r), alternate (a) we create upstream (up) and downstream (dn) kmers.
So the columns are all prefixed with one of "rup_", "rdn_", "aup_", "adn_".

In the alignment information per-query, there are three 'hit' counts:
- nhits : number of query hits
- dir_hits : direct strand hit count
- com_hits : compliment strand hit count

The rest of the alignment information is reported by average (avg), maximum (max), and minimum (min)

The suffixes are:
- q : mapping quality score of the hits
- ed : edit distance of the hits
- mat : number of matches
- mis : number of mismatches

For example, "aup_avg_q", is the alternate's upstream breakend kmer's average mapping quality score.

```
usage: grm [-h] -i INPUT -r REFERENCE [-R REGIONS] [-o OUTPUT] [-k KMERSIZE] [-m MIN_SIZE] [-t THREADS] [--debug]

Maps graph edge kmers with BWA to assess Graph Reference Mappability

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input VCF
  -r REFERENCE, --reference REFERENCE
                        BWA indexed reference
  -R REGIONS, --regions REGIONS
                        Bed file of regions to parse (None)
  -o OUTPUT, --output OUTPUT
                        Output dataframe (results.jl)
  -k KMERSIZE, --kmersize KMERSIZE
                        Size of kmer to map (50)
  -m MIN_SIZE, --min-size MIN_SIZE
                        Minimum size of variants to map (25)
  -t THREADS, --threads THREADS
                        Number of threads (1)
  --debug               Verbose logging
```

# truvari anno repmask

```
usage: repmask [-h] -i INPUT [-o OUTPUT] [-e EXECUTABLE] [-m MIN_LENGTH] [-M MAX_LENGTH] [-t THRESHOLD] [-p PARAMS] [-T THREADS]
               [--debug]

 Wrapper around RepeatMasker to annotate insertion sequences in a VCF

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        VCF to annotate (None)
  -o OUTPUT, --output OUTPUT
                        Output filename (/dev/stdout)
  -e EXECUTABLE, --executable EXECUTABLE
                        Path to RepeatMasker (RepeatMasker)
  -m MIN_LENGTH, --min-length MIN_LENGTH
                        Minimum size of entry to annotate (50)
  -M MAX_LENGTH, --max-length MAX_LENGTH
                        Maximum size of entry to annotate (50000)
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold for pct of allele covered (0.8)
  -p PARAMS, --params PARAMS
                        Default parameter string to send to RepeatMasker (-pa {threads} -qq -e hmmer -species human -lcambig
                        -nocut -div 50 -no_id -s {fasta})
  -T THREADS, --threads THREADS
                        Number of threads to use (1)
  --debug               Verbose logging
```

# truvari anno remap

Taking the Allele’s sequence, remap it to the reference and annotate based on the closest alignment.

![](https://github.com/spiralgenetics/truvari/blob/develop/imgs/remap_example.png)

```
usage: remap [-h] -r REFERENCE [-o OUTPUT] [-m MINLENGTH] [-t THRESHOLD] [-d DIST] [-H HITS] [--debug] [input]

Remap VCF'S alleles sequence to the reference to annotate REMAP

- novel : Allele has no hits in reference
- tandem : Allele's closest hit is within len(allele) bp of the SV's position
- interspersed : Allele's closest hit is not tandem
- partial : Allele only has partial hit(s) less than --threshold

Which alleles and alignments to consider can be altered with:
- --minlength : minimum SV length to considred (50)
- --dist : For deletion SVs, do not consider alignments that hit within Nbp of the SV's position
(a.k.a. alignments back to the source sequence) (10)
- --threshold : Minimum percent of allele's sequence used by alignment to be considered (.8)

positional arguments:
  input                 Input VCF (/dev/stdin)

options:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        BWA indexed reference
  -o OUTPUT, --output OUTPUT
                        Output VCF (/dev/stdout)
  -m MINLENGTH, --minlength MINLENGTH
                        Smallest length of allele to remap (50)
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold for pct of allele covered to consider hit (0.8)
  -d DIST, --dist DIST  Minimum distance an alignment must be from a DEL's position to be considered (10))
  -H HITS, --hits HITS  Report top hits as chr:start-end.pct (max 0)
  --debug               Verbose logging
```
# truvari anno hompct

```
usage: hompct [-h] -i INPUT [-o OUTPUT] [-b BUFFER] [-m MINANNO] [-M MAXGT] [-c MINCOUNT] [--debug]

Calcluate the the Hom / (Het + Hom) of variants in the region of SVs
Requires the VCF to contain SVs beside SNPs/Indels

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Compressed, indexed VCF to annotate
  -o OUTPUT, --output OUTPUT
                        Output filename (stdout)
  -b BUFFER, --buffer BUFFER
                        Number of base-pairs up/dn-stream to query (5000)
  -m MINANNO, --minanno MINANNO
                        Minimum size of event to annotate (50)
  -M MAXGT, --maxgt MAXGT
                        Largest event size to count for genotyping (1)
  -c MINCOUNT, --mincount MINCOUNT
                        Minimum number of genotyping events to report HOMPCT (0)
  --debug               Verbose logging
```

# truvari anno numneigh

```
usage: numneigh [-h] [-o OUTPUT] [-r REFDIST] [-s SIZEMIN] [--passonly] [--debug] [input]

For every call within size boundaries,
Add NumNeighbors info field of how many calls are within the distance
Add NeighId clustering field in the same chained neighborhood
For example,
::
    -- is a call, refdist is 2
         - - -   -    - -
    nn:  1 2 1   0    1 1
    id:  0 0 0   1    2 2

positional arguments:
  input                 VCF to annotate

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output vcf (stdout)
  -r REFDIST, --refdist REFDIST
                        Max reference location distance (1000)
  -s SIZEMIN, --sizemin SIZEMIN
                        Minimum variant size to consider for annotation (50)
  --passonly            Only count calls with FILTER == PASS
  --debug               Verbose logging
```

# truvari anno svinfo

Uses `truvari.entry_size` and  `truvari.entry_variant_type` on entries >= `args.minsize` to add 'SVLEN' and ‘SVTYPE’ annotations to a VCF’s INFO.

How SVLEN is determined:
- Starts by trying to use INFO/SVLEN
- If SVLEN is unavailable and ALT field is an SV (e.g. \<INS\>, \<DEL\>, etc), use abs(vcf.start - vcf.end). The INFO/END tag needs to be available, especially for INS.
- Otherwise, return the size difference of the sequence resolved call using abs(len(vcf.REF) - len(str(vcf.ALT[0])))

How SVTYPE is determined:
- Starts by trying to use INFO/SVTYPE
- If SVTYPE is unavailable, infer if entry is a insertion or deletion by looking at the REF/ALT sequence size differences
- If REF/ALT sequences are not available, try to parse the \<INS\>, \<DEL\>, etc from the ALT column.
- Otherwise, assume 'UNK'

```
usage: svinfo [-h] [-o OUTPUT] [-m MINSIZE] [input]

Adds SVTYPE and SVLEN INFO fields

positional arguments:
  input                 VCF to annotate (stdin)

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output filename (stdout)
  -m MINSIZE, --minsize MINSIZE
                        Minimum size of entry to annotate (50)
```

# truvari anno bpovl

After turning a tab-delimited annotation file into an IntervalTree, intersect the start/end and overlap of SVs.
The output is a light-weight pandas DataFrame saved with joblib. The columns in the output are:

- vcf_key : Variant key from `truvari.entry_to_key`
- intersection : Type of intersection between the SV and the annotation
    - start_bnd - SV's start breakpoints hits the annotation
    - end_bnd - SV's end breakpoint hits the annotation
    - overlaps - Annotation's start/end boundaries are completely within the SV
    - contains - Annotation's start/end boundaries span the entire SV
- anno_key : Annotation file's line index

The idea with this tool is to annotate variants against tab-delimited files, especially when there's a 1-to-N variant to annotations. This tool is useful when used in conjunction with `truvari vcf2df` and pandas DataFrames.

For example, if we have a VCF of SVs and a GTF of genes/gene features from Ensmbl. Any SV may intersect multiple features, which doesn't lend itself well to directly annotating the VCF's INFO. After using `bpovl`, we'll use Truvari to convert the SVs to a DataFrame.

```bash
truvari anno bpovl -i variants.vcf.gz -a genes.gtf.gz -o anno.jl -p gff
truvari vcf2df variants.vcf.gz variants.jl
```

We can then treat the files similar to a database and do queries and joins to report which variants intersect which annotations.

```python
import joblib
from gtfparse import read_gtf
variants = joblib.load("variants.jl")
genes = read_gtf("genes.gtf.gz")
annos = joblib.load("anno.jl")
to_check = annos.iloc[0]

print(to_check)
# vcf_key         chr20:958486-958487.A
# intersection                start_bnd
# anno_key                           11

print(variants.loc[to_check['vcf_key']])
# id                        None
# svtype                     INS
# ... etc

print(annos.loc[to_check['anno_key']])
# seqname                               chr20
# source                       ensembl_havana
# feature                                exon
# start                                958452
# ... etc
```

Similar to tabix, `bpovl` has presets for known file types like bed and gff. But any tab-delimited file with sequence/chromosome, start position, and end position can be parsed. Just set the "Annotation File Arguments" to the 0-based column indexes. For example, a bed file 
has arguments `-s 0 -b 1 -e 2  -c #`.

```
usage: bpovl [-h] -a ANNO -o OUTPUT [--sizemin SIZEMIN] [--spanmax SPANMAX] [-p {bed,gff}] [-c COMMENT] [-s SEQUENCE] [-b BEGIN]
             [-e END] [-1]
             [input]

Creates intersection of features in an annotation file with SVs' breakpoints and overlap

positional arguments:
  input                 VCF to annotate (stdin)

options:
  -h, --help            show this help message and exit
  -a ANNO, --anno ANNO  Tab-delimited annotation file
  -o OUTPUT, --output OUTPUT
                        Output joblib DataFrame
  --sizemin SIZEMIN     Minimum size of variant to annotate (50)
  --spanmax SPANMAX     Maximum span of SVs to annotate (50000)

Annotation File Arguments:
  -p {bed,gff}, --preset {bed,gff}
                        Annotation format. This option overwrites -s, -b, -e, -c and -1 (None)
  -c COMMENT, --comment COMMENT
                        Skip lines started with character. (#)
  -s SEQUENCE, --sequence SEQUENCE
                        Column of sequence/chromosome name. (0)
  -b BEGIN, --begin BEGIN
                        Column of start chromosomal position. (1)
  -e END, --end END     Column of end chromosomal position. (2)
  -1, --one-based       The position in the anno file is 1-based rather than 0-based. (False)
```
# truvari anno density
Partitions a `--genome` into `--windowsize` regions and count how many variants overlap. Annotate
regions with no variants as 'sparse' and with greater than or equal to (mean + `--threshold` * standard
deviation) number of variants as 'dense'. Outputs a joblib DataFrame with columns 
`chrom, start, end, count, anno`.

```
usage: density [-h] -g GENOME -o OUTPUT [-m MASK] [-w WINDOWSIZE] [-s STEPSIZE] [-t THRESHOLD] [input]

Identify 'dense' and 'sparse' variant windows of the genome

positional arguments:
  input                 Input VCF (/dev/stdin)

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        Genome bed file
  -o OUTPUT, --output OUTPUT
                        Output joblib DataFrame
  -m MASK, --mask MASK  Mask bed file
  -w WINDOWSIZE, --windowsize WINDOWSIZE
                        Window size (10000)
  -s STEPSIZE, --stepsize STEPSIZE
                        Window step size (10000)
  -t THRESHOLD, --threshold THRESHOLD
                        std for identifying 'dense' regions (3)
```

# truvari anno dpcnt

For multi-sample VCFs, it is often useful to have summarized depth (DP) information across samples per-variant. This adds a `INFO/DPCNT` with counts of how many samples have `FORMAT/DP` for each of the user-defined bins. Bins are incremented using `bisect` e.g. `pos = bisect.bisect(bins, dp); bins[pos] += 1;

```
usage: dpcnt [-h] [-b BINS] [--no-ad] [-p] [-o OUTPUT] [input]

Quick utility to count how many samples have >= Nx coverage per-variant

positional arguments:
  input                 VCF to annotate (stdin)

options:
  -h, --help            show this help message and exit
  -b BINS, --bins BINS  Coverage bins to bisect left the counts (0,5,10,15)
  --no-ad               Skip adding ADCNT bins
  -p, --present         Only count sites with present (non ./.) genotypes
  -o OUTPUT, --output OUTPUT
                        Output filename (stdout)
```

# truvari anno lcr

```
usage: lcr [-h] [-o OUTPUT] [input]

Annotate low complexity region entropy score for variants
Credit: https://jszym.com/blog/dna_protein_complexity/

positional arguments:
  input                 VCF to annotate (stdin)

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output filename (stdout)
```

# truvari anno grpaf

Add INFO tags of allele frequency annotations for groups of samples. For every group in `--labels` tab-delimited file, calculate the AF,MAF,ExcHet,HWE,MAC,AC for the samples in the group. Adds INFO tags with suffix of the group identifier (e.g. `AF_EAS`). `--strict` will hard fail if there are samples in the `--labels` not present in the vcf header.

```
usage: grpaf [-h] [-o OUTPUT] -l LABELS [-t TAGS] [--strict] [--debug] [input]

Add allele frequency annotations for subsets of samples

positional arguments:
  input                 VCF to annotate

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output filename (stdout)
  -l LABELS, --labels LABELS
                        Tab-delimited file of sample and group
  -t TAGS, --tags TAGS  Comma-separated list of tags to add from AF,MAF,ExcHet,HWE,MAC,AC,AN (all)
  --strict              Exit if sample listed in labels is not present in VCF (False)
  --debug               Verbose logging
```