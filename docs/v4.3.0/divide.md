Divide a VCF into independent shards.

Unfortunately, `pysam.VariantRecord` objects aren't pickle-able. This means that if we wanted to have Truvari leverage python's multiprocessing we'd need to make a custom VCF parser. However, the command `truvari divide` allows us to take an input VCF and divide it into multiple independent parts (or shards) which can the be processed over multiple processes. 

`truvari divide` works by parsing a VCF and splitting it into multiple, smaller sub-VCFs. If any variants are within `--buffer` base-pairs, they're output to the same sub-VCF. This allows variants in the same region which would need to be compared to one-another (see `--refdist`) to stay in the same sub-VCF. The `--min` parameter allows us to control the minimum number of variants per sub-VCF so that we don't make too many tiny VCFs. Once the sub-VCFs are created, we can process each independently through whatever truvari command.

For example, let's say we want to run `truvari collapse` on a very large VCF with many variants and many samples. First, we divide the VCF:

```bash
truvari divide big_input.vcf.gz sub_vcfs_directory/
```

Inside of `sub_vcfs_directory/` we'll have multiple VCFs, which we can process with a simple bash script

```bash
NJOBS=$(nproc) # use all processors by default
mkdir -p output_vcfs/
mkdir -p collap_vcfs/
mkdir -p logs/

for in_vcf in sub_vcfs_directory/*.vcf.gz
do
    # Setup file names
    base_name=$(basename $in_vcf)
    base_name=${base_name%.vcf.gz}
    output_name=output_vcfs/${base_name}.vcf
    collap_name=collap_vcfs/${base_name}.vcf
    log_name=logs/${base_name}.log
    # Run the command and send it to the background
    truvari collapse -i $in_vcf -o $output_name -c $collap_name -f reference.fa &> logs/${log_name}.log &
    # If too many jobs are running, wait for one to finish
    while [ $( jobs | wc -l ) -ge ${NJOBS} ]
    do
        sleep 5
    done
done
```

Obviously the logs and `while` loop are tricks for running on a single machine. If you have access to a cluster, I'm sure you can imagine how to create/submit the commands.

```
usage: divide [-h] [-b BUFFER] [-m MIN] [--no-compress] [-T THREADS] VCF DIR

Divide a VCF into independent parts

positional arguments:
  VCF                   VCF to split
  DIR                   Output directory to save parts

options:
  -h, --help            show this help message and exit
  -b BUFFER, --buffer BUFFER
                        Buffer to make mini-clusters (1000)
  -m MIN, --min MIN     Minimum number of entries per-vcf (100)
  --no-compress         Don't attempt to compress/index sub-VCFs
  -T THREADS, --threads THREADS
                        Number of threads for compressing/indexing sub-VCFs (1)
```