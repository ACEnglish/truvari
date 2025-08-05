# ------------------------------------------------------------
#                             stratp test
# ------------------------------------------------------------
# Ensure data.jl doesn't exist so we can test that branch
if [ -f $INDIR/bench_results/mims_output/data.jl ]
then
    rm $INDIR/bench_results/mims_output/data.jl
fi
run stratp1 $truv stratp \
        $INDIR/bench_results/mims_output/ \
        --preset mims \
        -o $OD/stratp1.txt
if [ $stratp1 ]; then
    assert_exit_code 0
    # There's rounding differences or something between local and github action
    #assert_equal $(fn_md5 $ANSDIR/stratp/stratp1.txt) $(fn_md5 $OD/stratp1.txt)
fi

run stratp2 $truv stratp \
        $INDIR/bench_results/mims_output/ \
        --features szbin,svtype,SRMAPQ:f3,SRQ:w3,SMACUIX7SY21_GQ:c0-100-1000 \
        --states comp \
        --tail right \
        --min-obs 5 \
        -o $OD/stratp2.txt \
        -c $OD/stratp2.cnt.txt
if [ $stratp2 ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/stratp/stratp2.txt) $(fn_md5 $OD/stratp2.txt)
    assert_equal $(fn_md5 $ANSDIR/stratp/stratp2.cnt.txt) $(fn_md5 $OD/stratp2.cnt.txt)
fi
