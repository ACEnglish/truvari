DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#For an explination of this, see https://github.com/spiralgenetics/truvari/wiki/Comparing-two-SV-programs
cmp_A=$1
cmp_B=$2

cd $cmp_A
paste <(grep -v "#" tp-base.vcf) <(grep -v "#" tp-call.vcf) > combined_tps.txt
cd ../$cmp_B
paste <(grep -v "#" tp-base.vcf) <(grep -v "#" tp-call.vcf) > combined_tps.txt
cd ../

(grep -v "#" $cmp_A/fn.vcf && grep -v "#" $cmp_B/fn.vcf) | cut -f3 | sort | uniq -c | grep "^ *1 " | cut -f2- -d1 > missed_names.txt

cat missed_names.txt | xargs -I {} grep -w {} $cmp_A/combined_tps.txt > missed_by_B.txt
cat missed_names.txt | xargs -I {} grep -w {} $cmp_B/combined_tps.txt > missed_by_A.txt

python $DIR/multi_venn.py $cmp_A/tp-base.vcf $cmp_B/tp-base.vcf

