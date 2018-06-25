import os
import sys
import json

from collections import defaultdict, Counter

import vcf
from biograph import cmd_exe

#TODO: document what this thing actually does
#TODO: pull cmd_exe over
#Bug fixes:
# Fixed genotype comparisons bug
# Spelling errors in documentation
# Added safety checks on ensuring input vcfs exist
# Added debug statements for WHY variants don't match. Should help trouble shoot why
#   events you expect to match don't and help users address the problem
# Added a script `benchmark_pipe.py` That will characterize snps/indels/svs in multiple
#   ways so that you can get a full characterization of performance (only relevant to
#   callers that report all variation regardless of size.)
# Added a link in README to truvari Wiki so it's more likely people will read that
#   extra information

# Separate SNPs/INDELs (just parse the tp/fns? or split and re-run...?)
# report those individual stats
# and I've already got giab-report made.

def do_cnt(bdir):
    cnt = {'tp': Counter(), 'fp': Counter(), 'fn': Counter()}
    for fn, label in [("tp.vcf.gz", "tp"), ("fp.vcf.gz", "fp"), ("fn.vcf.gz", "fn")]:
        v = vcf.Reader(filename=os.path.join(bdir, fn))
        for entry in v:
            cnt[label][entry.var_subtype] += 1
    return cnt


def summarize_cnt(cnt, out_file):
    with open(out_file, 'w') as fout:
        types = set(cnt['tp'].keys()).union(set(cnt['fp'].keys())).union(set(cnt['fn'].keys()))
        for cur_type in types:
            denom = float(cnt["tp"][cur_type]) + cnt["fn"][cur_type]
            if denom == 0:
                cnt[cur_type + "_recall"] = 0
            else:
                cnt[cur_type + "_recall"] = cnt["tp"][cur_type] / denom
            denom = float(cnt["tp"][cur_type]) + cnt["fp"][cur_type]
            if denom == 0:
                cnt[cur_type + "_precision"] = 0
            else:
                cnt[cur_type + "_precision"] = cnt["tp"][cur_type] / denom
        row_headers = []
        col_headers = None
        rows = []
        extra = []
        for key in cnt:
            if type(cnt[key]) is dict:
                row_headers.append(k)
                if col_headers is None:
                    col_headers = cnt[key].keys()
                c_row  = []
                for col in col_headers:
                    c_row.append(cnt[key][col])
                rows.append(c_row)
            else:
                extra.append("%s\t%d" % (key, cnt[key]))
        fout.write("\t" + "\t".join(col_headers) + '\n')
        for name, data in zip(row_headers, rows):
            fout.write("%s\t%s" % (name, "\t".join([str(x) for x in data])) + '\n')


def run():
    rtg_cmd = "rtg vcfeval -b {base_small} \
        -c {calls} \
        {high_conf} \
        {w_gt} \
        -t /share/datasets/HG001/hg19.sdf --all-records \
        -o results_small{hc}"
    truvari_cmd = "../truvari.py --base {base_large} \
        --comp {calls} \
        --giabreport \
        {w_gt} \
        {passonly} \
        {high_conf} \
        -o results_large{hc}"

    bdir = "/home/english/ajtrio/giab_calls"
    BASE_SMALL = os.path.join(
        bdir, "HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz")
    SMALL_HC = os.path.join(
        bdir, "HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed")
    BASE_LARGE = os.path.join(bdir, "HG002_SVs_Tier1_v0.6.vcf.gz")
    LARGE_HC = os.path.join(bdir, "HG002_SVs_Tier1_v0.6.bed")

    comp = sys.argv[1]
    
    print "HC small calls"
    r, o, e, t = cmd_exe(
        rtg_cmd.format(base_small=BASE_SMALL, high_conf="-e " + SMALL_HC, calls=comp, hc="_hc", w_gt="--squash-ploidy"))
    if r != 0:
        print r, o, e, t
    print o

    print "All small calls"
    r, o, e, c = cmd_exe(rtg_cmd.format(base_small=BASE_SMALL, high_conf="", calls=comp, hc="_all",
                                        w_gt="--squash-ploidy"))
    if r != 0:
        print r, o, e, t
    print o

    print "HC small calls w/GT"
    r, o, e, t = cmd_exe(
        rtg_cmd.format(base_small=BASE_SMALL, high_conf="-e " + SMALL_HC, calls=comp, hc="_hc_gt", w_gt=""))
    if r != 0:
        print r, o, e, c
    print o
    
    hc_cnt = do_cnt("results_small_hc")
    summarize_cnt(hc_cnt, "results_small_hc/by_type_cnt.txt")
    all_cnt = do_cnt("results_small_all")
    summarize_cnt(all_cnt, "results_small_all/by_type_cnt.txt")
    gt_cnt = do_cnt("results_small_hc_gt")
    summarize_cnt(gt_cnt, "results_small_hc_gt/by_type_cnt.txt")

    #checking SVs
    print "HC large calls"
    r, o, e, t = cmd_exe(truvari_cmd.format(base_large=BASE_LARGE, high_conf="--includebed " + LARGE_HC, calls=comp,
                                            hc="_hc", w_gt="", passonly="--passonly"))
    if r != 0:
        print r, o, e, t
    print o

    print "HC large calls w/ loose cmp"
    r, o, e, t = cmd_exe(truvari_cmd.format(base_large=BASE_LARGE, high_conf="--includebed " + LARGE_HC, calls=comp,
                                            hc="_hc_loose", w_gt="", passonly="--passonly --refdist 2000 --pctsim 0"))
    if r != 0:
        print r, o, e, t
    print o


    print "All large calls"
    r, o, e, c = cmd_exe(truvari_cmd.format(base_large=BASE_LARGE, high_conf="", calls=comp, hc="_all", w_gt="",
                passonly=""))
    if r != 0:
        print r, o, e, t
    print o
    
    print "HC large calls w/GT"
    r, o, e, t = cmd_exe(truvari_cmd.format(base_large=BASE_LARGE, high_conf="--includebed " + LARGE_HC, calls=comp,
                                            hc="_hc_gt", w_gt="--gtcomp", passonly=""))
    if r != 0:
        print r, o, e, t
    print o
    
if __name__ == '__main__':
    run()
