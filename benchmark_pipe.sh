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

