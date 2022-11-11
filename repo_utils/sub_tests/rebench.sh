# ------------------------------------------------------------
#                                 rebench
# ------------------------------------------------------------


export PATH=$INDIR/external/fake_mafft/:$PATH 

run test_rebench $truv bench -b $INDIR/rebench_data/hg002_base.vcf.gz \
                 -c $INDIR/rebench_data/hg002_comp.vcf.gz \
                 --includebed $INDIR/rebench_data/h1_hc_tr_hg002.bed \
                 -s 5 \
                 -o $OD/rebench_output
run test_rebench $truv rebench -u -f $INDIR/rebench_data/chr20.fa.gz $OD/rebench_output
if [ $test_rebench ]; then
    assert_exit_code 0
    run test_rebench_result
    assert_equal $(fn_md5 $INDIR/rebench_data/bench_output_small/rebench.counts.txt) \
                 $(fn_md5 $OD/rebench_output/rebench.counts.txt)
    assert_equal $(fn_md5 $INDIR/rebench_data/bench_output_small/rebench.summary.json) \
                 $(fn_md5 $OD/rebench_output/rebench.summary.json)
fi
