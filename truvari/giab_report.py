"""
When running against GIAB SVs, we can make reports
"""
# pylint: disable=no-member
import os
import logging
from collections import defaultdict, Counter

import pysam

def make_giabreport(args, stats_box):
    """
    Create summaries of the TPs/FNs
    """
    def make_entries(vcf_file):
        """
        Turn all vcf entries into a list of dicts for easy parsing
        """
        ret = []
        vcf_reader = pysam.VariantFile(vcf_file)
        for entry in vcf_reader:
            data = dict(entry.info)
            for sample in entry.samples:
                val = list(entry.samples[sample]["GT"]) + ['None', 'None']
                data[sample + "_GT"] = "{0}/{1}".format(*val)
            data["CHROM"] = entry.chrom
            data["POS"] = entry.pos
            data["ID"] = entry.id
            if entry.filter is not None:
                data["FILTER"] = ";".join(entry.filter.keys())
            data["QUAL"] = entry.qual
            ret.append(data)
        return ret

    def count_by(key, docs, out):
        """
        get a count for all the keys
        key is a dictionary of count and their order
        """
        main_key = next(iter(key.keys()))
        cnt = Counter()
        for x in docs:
            cnt[x[main_key]] += 1
        for k in key[main_key]:
            out.write("%s\t%s\n" % (k, cnt[k]))

    def twoxtable(key1, key2, docs, out):
        """
        Parse any set of docs and create a 2x2 table
        """
        main_key1 = next(iter(key1.keys()))
        main_key2 = next(iter(key2.keys()))
        cnt = defaultdict(Counter)
        for x in docs:
            cnt[x[main_key1]][x[main_key2]] += 1

        out.write(".\t" + "\t".join([str(i) for i in key1[main_key1]]) + '\n')
        for y in key2[main_key2]:
            o = [str(y)]
            for x in key1[main_key1]:
                o.append(str(cnt[x][y]))
            out.write("\t".join(o) + '\n')

    def bool_counter(keys, docs, out):
        """
        For bool valued INFO keys, write summaries of their counts
        """
        cnt = defaultdict(Counter)
        col_names = set()
        for x in docs:
            for y in keys:
                cnt[y][x[y]] += 1
                col_names.add(x[y])
        col_names = sorted(list(col_names))
        out.write(".\t%s\n" % "\t".join(col_names))
        row_keys = sorted(list(cnt.keys()))
        for row in row_keys:
            out.write(row)
            for col_val in col_names:
                out.write("\t%d" % cnt[row][col_val])
            out.write("\n")

    def collapse_techs(docs):
        """
        Make a new annotation about the presence of techs inplace
        called techs
        Illcalls
            PBcalls
        CGcalls
        TenXcalls
        """
        calls = ["Illcalls", "PBcalls", "CGcalls", "TenXcalls"]
        for d in docs:
            new_anno = []
            for i in calls:
                if d[i] > 0:
                    new_anno.append(i[:-5]) # .rstrip("calls"))
            d["techs"] = "+".join(new_anno)

    logging.info("Creating GIAB report")

    sum_out = open(os.path.join(args.output, "giab_report.txt"), 'w')

    tp_base = make_entries(os.path.join(args.output, "tp-base.vcf"))

    collapse_techs(tp_base)
    fn = make_entries(os.path.join(args.output, "fn.vcf"))
    collapse_techs(fn)

    size_keys = {"sizecat": ["50to99", "100to299", "300to999", "gt1000"]}
    svtype_keys = {"SVTYPE": ["DEL", "INS", "COMPLEX"]}
    tech_keys = {"techs": ["I+PB+CG+TenX", "I+PB+CG", "I+PB+TenX", "PB+CG+TenX",
                           "I+PB", "I+CG", "I+TenX", "PB+CG", "PB+TenX", "CG+TenX",
                           "I", "PB", "CG", "TenX"]}
    rep_keys = {"REPTYPE": ["SIMPLEDEL", "SIMPLEINS", "DUP", "SUBSDEL", "SUBSINS", "CONTRAC"]}
    tr_keys = ["TRall", "TRgt100", "TRgt10k", "segdup"]

    gt_keys_proband = {"HG002_GT": ["0/1", "./1", "1/1"]}
    gt_keys_father = {"HG003_GT": ["./.", "0/0", "0/1", "./1", "1/1"]}
    gt_keys_mother = {"HG004_GT": ["./.", "0/0", "0/1", "./1", "1/1"]}
    # OverallNumbers
    sum_out.write("TP\t%s\n" % (len(tp_base)))
    sum_out.write("FN\t%s\n\n" % (len(fn)))
    sum_out.write("TP_size\n")
    count_by(size_keys, tp_base, sum_out)
    sum_out.write("FN_size\n")
    count_by(size_keys, fn, sum_out)
    sum_out.write("\nTP_type\n")
    count_by(svtype_keys, tp_base, sum_out)
    sum_out.write("FN_type\n")
    count_by(svtype_keys, fn, sum_out)
    sum_out.write("\nTP_Type+Size\n")
    twoxtable(svtype_keys, size_keys, tp_base, sum_out)
    sum_out.write("FN_Type+Size\n")
    twoxtable(svtype_keys, size_keys, fn, sum_out)
    sum_out.write("\nTP_REPTYPE\n")
    count_by(rep_keys, tp_base, sum_out)
    sum_out.write("FN_REPTYPE\n")
    count_by(rep_keys, fn, sum_out)
    sum_out.write("\nTP_size+REPTYPE\n")
    twoxtable(size_keys, rep_keys, tp_base, sum_out)
    sum_out.write("FN_size+REPTYPE\n")
    twoxtable(size_keys, rep_keys, fn, sum_out)
    sum_out.write("\nTP_Tech\n")
    count_by(tech_keys, tp_base, sum_out)
    sum_out.write("FN_Tech\n")
    count_by(tech_keys, fn, sum_out)
    sum_out.write("\nTP_Size+Tech\n")
    twoxtable(size_keys, tech_keys, tp_base, sum_out)
    sum_out.write("FN_Size+Tech\n")
    twoxtable(size_keys, tech_keys, fn, sum_out)
    sum_out.write("\nTP_Type+Tech\n")
    twoxtable(svtype_keys, tech_keys, tp_base, sum_out)
    sum_out.write("FN_Type+Tech\n")
    twoxtable(svtype_keys, tech_keys, fn, sum_out)

    sum_out.write("\nPerformance\n")
    # Add output of the stats box
    for key in sorted(stats_box.keys()):
        sum_out.write("%s\t%s\n" % (key, str(stats_box[key])))

    sum_out.write("\nArgs\n")
    # Add output of the parameters
    argd = vars(args)
    for key in sorted(argd.keys()):
        sum_out.write("%s\t%s\n" % (key, str(argd[key])))
    #
    sum_out.write("\nTP_HG002GT\n")
    count_by(gt_keys_proband, tp_base, sum_out)
    sum_out.write("FN_HG002GT\n")
    count_by(gt_keys_proband, fn, sum_out)

    sum_out.write("\nTP_HG003.HG004GT\n")
    twoxtable(gt_keys_father, gt_keys_mother, tp_base, sum_out)
    sum_out.write("FN_HG003.HG004GT\n")
    twoxtable(gt_keys_father, gt_keys_mother, fn, sum_out)

    sum_out.write("\nTP TandemRepeat Anno\n")
    bool_counter(tr_keys, tp_base, sum_out)
    sum_out.write("FN TandemRepeat Anno\n")
    bool_counter(tr_keys, fn, sum_out)
    sum_out.close()
