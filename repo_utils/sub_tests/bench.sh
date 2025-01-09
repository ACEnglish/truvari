# ------------------------------------------------------------
#                                 bench
# ------------------------------------------------------------
bench() {
    # run and test truvari bench
    f1=$1
    f2=$2
    k=$3
    rm -rf $OD/bench${k}
    $truv bench -b $INDIR/variants/input${f1}.vcf.gz \
                -c $INDIR/variants/input${f2}.vcf.gz \
                -f $INDIR/references/reference.fa \
                --dup-to-ins \
                -o $OD/bench${k}/ ${4}
}

bench_assert() {
    k=$1
    assert_exit_code 0
    for i in $ANSDIR/bench/bench${k}/*.vcf.gz
    do
        bname=$(basename $i | sed 's/[\.|\-]/_/g')
        result=$OD/bench${k}/$(basename $i)
        assert_equal $(fn_md5 $i) $(fn_md5 $result)
    done
}

run bench_12 bench 1 2 12
if [ $bench_12 ]; then
    bench_assert  12
fi

run bench_13 bench 1 3 13 "--write-resolved"
if [ $bench_13 ]; then
    bench_assert 13
fi

run bench_23 bench 2 3 23 "--pick multi"
if [ $bench_23 ]; then
    bench_assert 23
fi

# Testing --includebed
run bench_13_includebed bench 1 3 13_includebed "--includebed $INDIR/beds/include.bed"
if [ $bench_13_includebed ]; then
    bench_assert 13_includebed
fi

# Testing --extend
run bench_13_extend bench 1 3 13_extend "--includebed $INDIR/beds/include.bed --extend 500"
if [ $bench_13_extend ]; then
    bench_assert 13_extend
fi


# --unroll
run bench_unroll $truv bench -b $INDIR/variants/real_small_base.vcf.gz \
                                  -c $INDIR/variants/real_small_comp.vcf.gz \
                                  -o $OD/bench_unroll/
if [ $bench_unroll ]; then
    bench_assert _unroll
fi

# --pick allele count
run bench_12_gtcomp bench 1 2 12_gtcomp "--pick ac"
if [ $bench_12_gtcomp ]; then
    bench_assert 12_gtcomp
fi

# --gtcomp edgecase
run bench_gtcomp_edgecase1 $truv bench -b $INDIR/variants/gtcomp_problem1_base.vcf.gz \
                                            -c $INDIR/variants/gtcomp_problem1_comp.vcf.gz \
                                            --pick ac \
                                            -o $OD/bench_gtcomp_edgecase1/
if [ $bench_gtcomp_edgecase1 ]; then
    bench_assert _gtcomp_edgecase1
fi
run bench_badparams $truv bench -b nofile.vcf -c nofile.aga -f notref.fa -o $OD --refdist 0 --extend
if [ $bench_badparams ]; then
    assert_exit_code 100
fi

run bench_starallele $truv bench -b $INDIR/variants/star.base.vcf.gz \
                                      -c $INDIR/variants/star.comp.vcf.gz \
                                      -s 0 -o $OD/bench_starallele/
if [ $bench_starallele ]; then
    bench_assert _starallele
fi

run bench_bnd $truv bench -b $INDIR/variants/bnd.base.vcf.gz \
                               -c $INDIR/variants/bnd.comp.vcf.gz \
                               -p 0 -o $OD/bench_bnd/ --no-decompose
if [ $bench_bnd ]; then
    bench_assert _bnd
fi

run bench_bnd_decomp $truv bench -b $INDIR/variants/bnd.base.vcf.gz \
                                 -c $INDIR/variants/bnd.comp2.vcf.gz \
                                 --sizemax 1000000000 -p 0 --pick multi \
                                 -o $OD/bench_bnd_decomp/
if [ $bench_bnd_decomp ]; then
    bench_assert _bnd_decomp/
fi
