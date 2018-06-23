DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "INCOMPLETE - do not use"
exit 1
#For an explination of this, see https://github.com/spiralgenetics/truvari/wiki/..?

base=$1
comp=$2
outdir=$3

#Loose Match Parameters
lm_param="-b $base -c $comp --pctsim 0 --pctsize 0.5 --refdist 1000"

#Allele Match Parameters
am_param="-b $base -c $comp --passonly"

#Genotype Match Parameters
gm_param="-b $base -c $comp --passonly --genotype"

#Run the individual programs
mkdir -f $outdir
python $DIR/truvari.py $lm_param -o $outdir/lm_out
python $DIR/truvari.py $lm_param -o $outdir/am_out
python $DIR/truvari.py $lm_param -o $outdir/gm_out


python $DIR/benchmark_consolidate.py $outdir

# collect gm_out/tp* and annotate while keeping track of variants used
# collect am_out/tp* and annotate those not already matched whlie keeping track of new variants used
# collect lm_out/tp* and annotate those not already matched while keeping track of new variants used
# for every $base call - if not used - annotate it as a failure (with reason?)
# for every $comp call - if not used - annotate it as a failure (with reason?)
# May also need to put on stratification annotation? - I might just make a tool to add them
# to a vcf. Then you can separate accordingly.
# So we'll add the stratification annotations to base/comp and benchmark_consolidate
#   will look for the results as part of the counting.
# Is Location in the TableOutput below count one of the stratification?
# https://github.com/ga4gh/benchmarking-tools/tree/master/resources/stratification-bed-files

# Example Table output?
#                        Metric 1    Metric 2    ...
#Location/count type 1   < value >   < value >   ...
#Location/count type 2   < value >   < value >   ...
#                  ...   < value >   < value >   ...


# The output VCF file is similar to the intermediate VCF file. We add one additional variant classification: in addition
# to FP/FN/TP, each variant call can also be assigned the status UNK for unknown / outside the regions which the
# truthset covers.

# A simple definition for UNK variants is as follows: We call any variant unknown if it BT == FP and BK == miss and if
# the variant is outside the confident call regions of the truth set. If the truthset does not give confident call
# regions, no UNK variants are output.

# Annotation 1
##FORMAT=<ID=BK,Number=1,Type=String,Description="Sub-type for decision (match/mismatch type)">
#The value of BK specifies the class of match for each variant record:

#.: missing = no match at any level tested by the comparison tool
#lm: loose match = the truth/query variant is nearby a variant in the query/truth -- if the tool outputs such match types, it should annotate the VCF header with the definition of loose matches (e.g. match within a fixed window, or within the same superlocus).
#am: almatch = the variant forms (part of) an allele match (independent of representation, i.e. one-sided haplotype match)
#gm: gtmatch = diploid haplotypes (and genotypes) were resolved to be the same (independent of representation)
##FORMAT=<ID=BD,Number=1,Type=String,Description="Decision for call (TP/FP/FN/N)">

