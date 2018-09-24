import os
import sys
import time
import json
import signal
import datetime
import argparse
import subprocess

from collections import defaultdict, Counter

import vcf

USAGE = """Run many stratifications of SNP/INDEL/SV benchmarking. 

If small variants (-a and -A) are not provided, rtg will not be run"""

class Alarm(Exception):
    """ Alarm Class for command timeouts """
    def __init__(self):
        pass


def alarm_handler(signum, frame=None):  # pylint: disable=unused-argument
    """ Alarm handler for command timeouts """
    raise Alarm


def cmd_exe(cmd, timeout=-1):
    """
    Executes a command through the shell.
    timeout in minutes! so 1440 mean is 24 hours.
    -1 means never
    returns (ret_code, stdout, stderr, datetime)
    where ret_code is the exit code for the command executed
    stdout/err is the Standard Output Error from the command
    and datetime is a datetime object of the execution time
    """
    t_start = time.time()
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, close_fds=True,
                            preexec_fn=os.setsid)
    signal.signal(signal.SIGALRM, alarm_handler)
    if timeout > 0:
        signal.alarm(int(timeout * 60))
    try:
        stdoutVal, stderrVal = proc.communicate()
        signal.alarm(0)  # reset the alarm
    except Alarm:
        logging.error(("Command was taking too long. "
                       "Automatic Timeout Initiated after %d"), timeout)
        os.killpg(proc.pid, signal.SIGTERM)
        proc.kill()
        return 214, None, None, datetime.timedelta(seconds=int(time.time() - t_start))
    t_end = time.time()

    stdoutVal = bytes.decode(stdoutVal)
    retCode = proc.returncode
    return retCode, stdoutVal, stderrVal, datetime.timedelta(seconds=int(t_end - t_start))


def do_cnt(bdir):
    cnt = {'tp': Counter(), 'fp': Counter(), 'fn': Counter()}
    for fn, label in [("tp.vcf.gz", "tp"), ("fp.vcf.gz", "fp"), ("fn.vcf.gz", "fn")]:
        v = vcf.Reader(filename=os.path.join(bdir, fn))
        for entry in v:
            cnt[label][str(entry.var_subtype)] += 1
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
            if type(cnt[key]) is Counter:
                row_headers.append(key)
                if col_headers is None:
                    col_headers = cnt[key].keys()
                c_row  = []
                for col in col_headers:
                    c_row.append(cnt[key][col])
                rows.append(c_row)
            else:
                extra.append("%s\t%s" % (key, str(cnt[key])))
        fout.write("\t" + "\t".join(col_headers) + '\n')
        for name, data in zip(row_headers, rows):
            fout.write("%s\t%s" % (name, "\t".join([str(x) for x in data])) + '\n')

def parse_args(args):
    """ Make pretty arguments """
    parser = argparse.ArgumentParser(description=USAGE)
    parser.add_argument("-a", "--base-small", required=False, default=None,
                        help="Small benchmarking base variants")
    parser.add_argument("-A", "--base-small-hc-region", required=False, default=None,
                        help="High Confidence Regions for small variants")
    parser.add_argument("-b", "--base-large", required=True,
                        help="Large benchmarking base variants")
    parser.add_argument("-B", "--base-large-hc-region", required=True,
                        help="High Confidence Regions for large variants")
    parser.add_argument("-c", "--calls", required=True,
                        help="Comparison variants")
    parser.add_argument("-t", "--rtg_ref", required=False,
                        help="RTG indexed reference")
    parser.add_argument("-r", "--reference", required=True,
                        help="samtools faidx indexed reference")
    parser.add_argument("-s", "--simple", action="store_true", default=False,
                        help="Only create the small/large hc results")
    args = parser.parse_args()
    return args

def run(args):
    args = parse_args(args)
    rtg_cmd = "rtg vcfeval -b {base_small} \
        -c {calls} \
        {high_conf} \
        {w_gt} \
        -t %s --all-records \
        -o results_small{hc}" % args.rtg_ref
    truvari_cmd = "truvari.py --base {base_large} \
        --comp {calls} \
        --reference %s \
        --giabreport \
        {w_gt} \
        {passonly} \
        {high_conf} \
        -o results_large{hc}" % args.reference

    BASE_SMALL = args.base_small
    SMALL_HC = args.base_small_hc_region
    BASE_LARGE = args.base_large
    LARGE_HC = args.base_large_hc_region
    comp = args.calls
    
    if args.rtg_ref is not None and BASE_SMALL is not None and SMALL_HC is not None:
        m_cmd = rtg_cmd.format(base_small=BASE_SMALL, high_conf="-e " + SMALL_HC, calls=comp, hc="_hc", w_gt="--squash-ploidy --sample HG002")
        print "HC small calls\n%s" % m_cmd
        r, o, e, t = cmd_exe(m_cmd)
            
        if r != 0:
            print r, o, e, t
        print o
        hc_cnt = do_cnt("results_small_hc")
        summarize_cnt(hc_cnt, "results_small_hc/by_type_cnt.txt")       

        if not args.simple:
            print "All small calls"
            r, o, e, c = cmd_exe(rtg_cmd.format(base_small=BASE_SMALL, high_conf="", calls=comp, hc="_all",
                                            w_gt="--squash-ploidy"))
            if r != 0:
                print r, o, e, t
            print o

            all_cnt = do_cnt("results_small_all")
            summarize_cnt(all_cnt, "results_small_all/by_type_cnt.txt")

            print "HC small calls w/GT"
            r, o, e, t = cmd_exe(
                rtg_cmd.format(base_small=BASE_SMALL, high_conf="-e " + SMALL_HC, calls=comp, hc="_hc_gt", w_gt=""))
            if r != 0:
                print r, o, e, c
            print o

            gt_cnt = do_cnt("results_small_hc_gt")
            summarize_cnt(gt_cnt, "results_small_hc_gt/by_type_cnt.txt")
    else:
        print "Skipping analysis of small events. Provide -a and -A and -t to run rtg"
    
    #checking SVs
    m_cmd = truvari_cmd.format(base_large=BASE_LARGE, high_conf="--includebed " + LARGE_HC, calls=comp, 
                                hc="_hc", w_gt="", passonly="--passonly")

    print "HC large calls\n%s" % m_cmd
    r, o, e, t = cmd_exe(m_cmd)
    if r != 0:
        print r, o, e, t
    print o
    
    if not args.simple:
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
                                            hc="_hc_gt", w_gt="--gtcomp", passonly="--passonly"))
        if r != 0:
            print r, o, e, t
        print o
    
if __name__ == '__main__':
    run(sys.argv)
