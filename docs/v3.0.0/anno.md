
Truvari annotations:
        gcpct, gtcnt, trf, grm, repmask, remap, hompct, numneigh, svinfo

# truvari anno gcpct

This will add an INFO tag `GCPCT` to each element in a VCF of the GC percent of the call's sequence.

For deletions, this is the GC percent of the reference range of the call. For insertions, the ALT sequence is analyzed.
```
usage: gcpct [-h] [-i INPUT] [-o OUTPUT] -r REFERENCE

Annotates GC content of SVs

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        VCF to annotate (stdin)
  -o OUTPUT, --output OUTPUT
                        Output filename (stdout)
  -r REFERENCE, --reference REFERENCE
                        Reference fasta
```

# truvari anno gtcnt
This will add an INFO tag `GTCNT` to each element in a VCF with the count of genotypes found across all samples. The value is a list of Counts of genotypes for the allele across all samples (UNK, REF, HET, HOM). This is most useful for pVCFs.

```
usage: gtcnt [-h] [-i INPUT] [-o OUTPUT]

Annotates GTCounts of alleles

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        VCF to annotate (stdin)
  -o OUTPUT, --output OUTPUT
                        Output filename (stdout)
```

# truvari anno trf
Adds a tandem-repeat annotation to each entry in your VCF.

This requires the [trf executable](http://tandem.bu.edu/trf/trf.download.html), which can either be in your environment `$PATH` or can be directly pointed to through the command line parameters. This also uses SimpleRepeat tracks over the reference genome (e.g. [`resources/`](https://github.com/spiralgenetics/truvari/tree/develop/resources)).

First, every SV is parsed and intersected with the SimpleRepeat track. If the SV hits a SimpleRepeat region, the SV will be annotated. Furthermore, we build the haplotype sequence of the alternate allele and the SimpleRepeat region and feed that sequence into tandem repeat finder. We then intersect the repeats TRF found with the repeats found in the SimpleRepeat track intersection. If the same repeat is found in both, we annotate with the TRF information as well as reporting the difference in the number of copies in alternate allele's haplotype sequence via TRFDiff.

INFO fields added
* TRF - Entry hits a simple repeat region
* TRFDiff - Simple repeat copy difference
* TRFperiod - Period size of the repeat
* TRFcopies - Number of copies aligned with the consensus pattern.
* TRFscore - TRF Alignment score
* TRFentropy - TRF Entropy measure
* TRFrepeat - TRF Repeat found on entry

```
usage: trf [-h] [-i INPUT] [-o OUTPUT] [-e EXECUTABLE] [-T TRF_PARAMS]
           [-s SIMPLE_REPEATS] [-f REFERENCE] [-m MIN_LENGTH] [-M MAX_LENGTH]
           [-t THREADS] [-C CHUNK_SIZE] [--debug]

Intersect vcf with reference simple repeats and report
how the an alternate allele affects the copies using TRF

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        VCF to annotate (stdin)
  -o OUTPUT, --output OUTPUT
                        Output filename (stdout)
  -e EXECUTABLE, --executable EXECUTABLE
                        Path to tandem repeat finder (trf409.linux64)
  -T TRF_PARAMS, --trf-params TRF_PARAMS
                        Default parameters to send to trf (3 7 7 80 5 40 500
                        -h -ngs)
  -s SIMPLE_REPEATS, --simple-repeats SIMPLE_REPEATS
                        Simple repeats bed
  -f REFERENCE, --reference REFERENCE
                        Reference fasta file
  -m MIN_LENGTH, --min-length MIN_LENGTH
                        Minimum size of entry to annotate (50)
  -M MAX_LENGTH, --max-length MAX_LENGTH
                        Maximum size of sequence to run through trf (10000)
  -t THREADS, --threads THREADS
                        Number of threads to use (48)
  -C CHUNK_SIZE, --chunk-size CHUNK_SIZE
                        Size (in mbs) of reference chunks for parallelization
                        (1)
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
usage: grm [-h] -i INPUT -r REFERENCE [-o OUTPUT] [-k KMERSIZE] [-m MIN_SIZE]
           [-t THREADS] [--debug]

Maps graph edge kmers with BWA to assess Graph Reference Mappability

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input VCF
  -r REFERENCE, --reference REFERENCE
                        BWA indexed reference
  -o OUTPUT, --output OUTPUT
                        Output dataframe (results.jl)
  -k KMERSIZE, --kmersize KMERSIZE
                        Size of kmer to map (50)
  -m MIN_SIZE, --min-size MIN_SIZE
                        Minimum size of variants to map (25)
  -t THREADS, --threads THREADS
                        Number of threads (48)
  --debug               Verbose logging
```

# truvari anno repmask

```
usage: repmask [-h] -i INPUT [-o OUTPUT] [-e EXECUTABLE] [-m MIN_LENGTH]
               [-M MAX_LENGTH] [-t THRESHOLD] [-p PARAMS] [-T THREADS]
               [--debug]

 Wrapper around RepeatMasker to annotate insertion sequences in a VCF

optional arguments:
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
                        Default parameter string to send to RepeatMasker (-pa
                        {threads} -e hmmer -species human -lcambig -nocut -div
                        50 -no_id -s {fasta})
  -T THREADS, --threads THREADS
                        Number of threads to use (48)
  --debug               Verbose logging
```

# truvari anno remap

Taking the Allele’s sequence, remap it to the reference and annotate based on the closest alignment.

![](https://github.com/spiralgenetics/truvari/blob/develop/imgs/remap_example.png)

```
usage: remap [-h] [-i INPUT] -r REFERENCE [-o OUTPUT] [-m MINLENGTH]
             [-t THRESHOLD] [-d DIST] [--debug]

Remap VCF'S alleles sequence to the reference to annotate REMAP

novel - Allele has no hits in reference
tandem - Allele's closest hit is within len(allele) bp of the SV's position
interspersed - Allele's closest hit is not tandem
partial - Allele only has partial hit(s) less than --threshold

Which alleles and alignments to consider can be altered with:
--minlength - minimum SV length to considred (50)
--dist - For deletion SVs, do not consider alignments that hit within Nbp of the SV's position
        (a.k.a. alignments back to the source sequence) (10)
--threshold - Minimum percent of allele's sequence used by alignment to be considered (.8)

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input VCF (/dev/stdin)
  -r REFERENCE, --reference REFERENCE
                        BWA indexed reference
  -o OUTPUT, --output OUTPUT
                        Output VCF (/dev/stdout)
  -m MINLENGTH, --minlength MINLENGTH
                        Smallest length of allele to remap (50)
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold for pct of allele covered to consider hit
                        (0.8)
  -d DIST, --dist DIST  Minimum distance an alignment must be from a DEL's
                        position to be considered (10))
  --debug               Verbose logging
```
# truvari anno hompct

```
usage: hompct [-h] -i INPUT [-o OUTPUT] [-b BUFFER] [-m MINANNO] [-M MAXGT]
              [--debug]

Calcluate the the Hom / (Het + Hom) of variants in the region of SVs
Requires the VCF to contain SVs beside SNPs/Indels

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Compressed, indexed VCF to annotate
  -o OUTPUT, --output OUTPUT
                        Output filename (stdout)
  -b BUFFER, --buffer BUFFER
                        Number of base-pairs up/dn-stream to query (5000)
  -m MINANNO, --minanno MINANNO
                        Minimum size of event to annotate
  -M MAXGT, --maxgt MAXGT
                        Largest event size to count for genotyping (1)
  --debug               Verbose logging
```

# truvari anno numneigh

```
usage: numneigh [-h] [-i INPUT] [-o OUTPUT] [-r REFDIST] [-s SIZEMIN]
                [--passonly] [--debug]

For every call within size boundaries,
Add NumNeighbors info field of how many calls are within the distance
Add NeighId clustering field in the same chained neighborhood
For example,
-- is a call, refdist is 2
     - - -   -    - -
nn:  1 2 1   0    1 1
id:  0 0 0   1    2 2

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        VCF to annotate
  -o OUTPUT, --output OUTPUT
                        Output vcf (stdout)
  -r REFDIST, --refdist REFDIST
                        Max reference location distance (1000)
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
usage: svinfo [-h] [-i INPUT] [-o OUTPUT] [-m MINSIZE]

Adds SVTYPE and SVLEN INFO fields

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        VCF to annotate (stdin)
  -o OUTPUT, --output OUTPUT
                        Output filename (stdout)
  -m MINSIZE, --minsize MINSIZE
                        Minimum size of entry to annotate (50)
  -s SIZEMIN, --sizemin SIZEMIN
                        Minimum variant size to consider for annotation (50)
  --passonly            Only count calls with FILTER == PASS
  --debug               Verbose logging
```