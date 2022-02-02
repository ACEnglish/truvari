"""
When running against GIAB SVs, we can make reports
"""
import os

import joblib
import pandas as pd
import truvari


def make_tech(data):
    """
    Add tech and tech_str columns to the data frame
    """
    flag = {1: "PB", 2: "Ill", 4: "TenX", 8: "CG"}
    data['tech'] = 0
    for val, col in enumerate(["PBcalls", "Illcalls", "TenXcalls", "CGcalls"]):
        data['tech'] += (data[col] != 0).astype(int) * (val + 1) ** 2

    def make_tech_str(x):
        ret = []
        for key, value in flag.items():
            if x & key:
                ret.append(value)
        if ret:
            return "+".join(ret)
        return None
    data["tech_str"] = data['tech'].apply(make_tech_str)


def make_giabreport(args, stats_box):
    """
    Create summaries of the TPs/FNs
    """
    out_df = os.path.join(args.output, "giab.jl")
    truvari.vcf2df.vcf2df_main(["-i", "-f", "-b", args.output, out_df])
    data = joblib.load(out_df)
    sz_type = pd.CategoricalDtype(
        categories=["50to99", "100to299", "300to999", "gt1000"], ordered=True)
    data["sizecat"] = data["sizecat"].astype(sz_type)
    make_tech(data)
    with open(os.path.join(args.output, "giab_report.txt"), 'w') as sum_out:
        sum_out.write("# Performance\n")
        sum_out.write(f"recall: {stats_box['recall']}\n")
        sum_out.write(f"precision: {stats_box['precision']}\n")
        sum_out.write(f"f1: {stats_box['f1']}\n")
        sum_out.write(f"gt_concordance: {stats_box['gt_concordance']}\n")
        sum_out.write(data["state"].value_counts().to_csv(sep="\t") + "\n")
        sum_out.write("\n")

        base_calls = data["state"].isin(["tpbase", "fn"])
        sum_out.write("# State by SizeCat and SVTYPE\n")
        view = data[base_calls].groupby(
            ["sizecat", "SVTYPE", "state"]).size().unstack(fill_value=0)
        view["pct"] = view["tpbase"] / view.sum(axis=1)
        sum_out.write(view.to_csv(sep="\t") + '\n')
        sum_out.write("\n")

        sum_out.write("# State by REPTYPE\n")
        view = data[base_calls].groupby(
            ["REPTYPE", "state"]).size().unstack(fill_value=0)
        view["pct"] = view["tpbase"] / (view["fn"] + view["tpbase"])
        view.sort_index(inplace=True)
        sum_out.write(view.to_csv(sep="\t") + '\n')
        sum_out.write("\n")

        sum_out.write("# State by Tech\n")
        view = data[base_calls].groupby(
            ["tech_str", "state"]).size().unstack(fill_value=0).sort_index()
        view["pct"] = view["tpbase"] / view.sum(axis=1)
        sum_out.write(view.to_csv(sep="\t") + '\n')
        sum_out.write("\n")

        sum_out.write("# State by Genotype\n")
        data["HG002_GT"] = data["GT"].apply(lambda x: truvari.get_gt(x).name)
        view = data.groupby(["state", "HG002_GT"]).size().unstack(fill_value=0)
        view.loc["total"] = view.sum(axis=0)
        view["TOT"] = view.sum(axis=1)
        view.loc["recall"] = view.loc["tpbase"] / \
            (view.loc["tpbase"] + view.loc["fn"])
        view.loc["precision"] = view.loc["tp"] / \
            (view.loc["tp"] + view.loc["fp"])
        view = view.round(decimals=2)
        sum_out.write(view.to_csv(sep="\t") + '\n')
        sum_out.write("\n")

        sum_out.write("# Args\n")
        argd = vars(args)
        for key in sorted(argd.keys()):
            sum_out.write(f"{key}\t{str(argd[key])}\n")
    joblib.dump(data, out_df)
