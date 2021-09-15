
Truvari annotations:
        gcpct, gtcnt, trf, grm, repmask, remap, hompct, numneigh

# truvari anno gcpct

This will add an INFO tag `GCPCT` to each element in a VCF of the GC percent of the call's sequence.

For deletions, this is the GC percent of the reference range of the call. For insertions, the ALT sequence is analyzed.
```
usage: truvari anno gcpct [-h] [-i INPUT] [-o OUTPUT] -r REFERENCE

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
usage: truvari anno gtcnt [-h] [-i INPUT] [-o OUTPUT]

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

This requires [trf executable](http://tandem.bu.edu/trf/trf.download.html), which can either be in your environment `$PATH` or can be directly pointed to through the command line parameters. Every entry will have the deleted or inserted sequence run through TRF, and the returned fields will be added to the entry's INFO fields.

Use the `--ref-bed` with a reference based truvari annotation track  (`resources/*simpleRepeat.truvari.txt.gz`). Any intersection to a known tandem repeat in the reference will also be added to the entry's INFO field. 

An additional INFO field `TRF_diffs` is populated where any tandem repeat motifs shared between the reference annotation and the TRF annotation over the entry's sequence are the same, we report the Alternate allele copy number difference from its reference SREP_repeat copy number. These numbers are reported 1-to-1 in the same order of `TRF_repeats`.

```
usage: truvari anno trf [-h] [-i INPUT] [-o OUTPUT] [-e EXECUTABLE] [-f] [-m MIN_LENGTH]
           [-b BESTN] [-t TRF_PARAMS] [-R REF_BED] [-r REF] [-l] [--debug]

 Wrapper around TRF to annotate a VCF 

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        VCF to annotate (stdin)
  -o OUTPUT, --output OUTPUT
                        Output filename (stdout)
  -e EXECUTABLE, --executable EXECUTABLE
                        Path to tandem repeat finder (trf409.linux64)
  -f, --full            Write full trf output to entries
  -m MIN_LENGTH, --min-length MIN_LENGTH
                        Minimum size of entry to annotate (50)
  -b BESTN, --bestn BESTN
                        Top trf hits to report (5)
  -t TRF_PARAMS, --trf-params TRF_PARAMS
                        Default parameters to send to tr (2 7 7 80 10 50 500
                        -m -f -h -d -ngs)
  -R REF_BED, --ref-bed REF_BED
                        Reference bed of tandem repeat regions
  -r REF, --ref REF     Reference fasta file (use with -R)
  -l, --layered         (non-functional)
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
usage: grm [-h] -i INPUT -r REFERENCE [-o OUTPUT] [-k KMERSIZE] [-t THREADS] [--debug]

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
  -t THREADS, --threads THREADS
                        Number of threads (4)
  --debug               Verbose logging
```

# truvari anno repmask

```
usage: repmask [-h] -i INPUT [-o OUTPUT] [-e EXECUTABLE] [-m MIN_LENGTH] [-M MAX_LENGTH] [-t THRESHOLD] [-p PARAMS] [-T THREADS] [--debug]

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
                        Default parameter string to send to RepeatMasker (-pa {threads} -e hmmer -species human -lcambig -nocut -div 50
                        -no_id -s {fasta})
  -T THREADS, --threads THREADS
                        Number of threads to use (4)
  --debug               Verbose logging
```

# truvari anno remap
Note that 'novel' in the context of Deletions means that the deleted sequence has no *other* hits in reference besides the source location.
```
usage: remap [-h] [-i INPUT] -r REFERENCE [-o OUTPUT] [-m MINLENGTH]
             [-t THRESHOLD] [-d DIST] [--debug]

Remap VCF'S alleles sequence to the reference to annotate REMAP

novel - Allele has no hits in reference
tandem - Allele's closest hit is within len(allele) bp of the SV's position
interspersed - Allele's closest hit is not tandem

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
                        Threshold for pct of allele covered (0.8)
  -d DIST, --dist DIST  Minimum distance an alignment must be from a DEL's
                        position to be considered (10))
  --debug               Verbose logging
```
# truvari anno hompct

```
usage: hompct [-h] -i INPUT [-o OUTPUT] [-b BUFFER] [-m MINANNO] [-M MAXGT] [--debug]

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
usage: numneigh [-h] -i INPUT [-o OUTPUT] [-r REFDIST] [-s SIZEMIN] [--passonly] [--debug]

For every call within size boundaries,
Add NumNeighbors info field of how many calls over within size boundary
are in the neighborhood

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        VCF to annotate
  -o OUTPUT, --output OUTPUT
                        Output vcf (stdout)
  -r REFDIST, --refdist REFDIST
                        Max reference location distance (1000)
  -s SIZEMIN, --sizemin SIZEMIN
                        Minimum variant size to consider for annotation (50)
  --passonly            Only count calls with FILTER == PASS
  --debug               Verbose logging
```