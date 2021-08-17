"""
Takes a vcf and creates a data frame. Can parse a bench output directory
"""
import os
import sys
import glob
import logging
import argparse

import pysam
import joblib
import pandas as pd
import truvari

SZBINTYPE = pd.CategoricalDtype(categories=truvari.SZBINS, ordered=True)
SVTYTYPE = pd.CategoricalDtype(categories=[_.name for _ in truvari.SV], ordered=True)

def get_files_from_truvdir(directory):
    """
    Finds one of each of the files in the directory we're expection
    Raises FileNotFoundError if they're not there or there are too many
    return dictionary of file identifier and path (e.g. "tp-base": "dir/tp-base.vcf.gz")
    """
    ret = {}
    pats = [("tpbase", "tp-base.vcf"), ("tp", "tp-call.vcf"),
            ("fp", "fp.vcf"), ("fn", "fn.vcf")]
    has_error = False
    for key, name in pats:
        fname = os.path.join(directory, name)
        ret[key] = glob.glob(fname)
        ret[key].extend(glob.glob(fname + '.gz'))
        if len(ret[key]) != 1:
            logging.error("Expected 1 file for %s[.gz] found %d", fname, len(ret[key]))
            has_error = True
    #summary = glob.glob(os.path.join(directory, "summary.txt"))
    #ret["summary"]
    if has_error:
        raise FileNotFoundError(f"Couldn't parse truvari directory {directory}. See above errors")
    return ret

def vcf_to_df(fn, with_info=True, with_fmt=True, sample=0):
    """
    Parse a vcf file and turn it into a dataframe.
    Tries its best to pull info/format tags from the sample information
    For Formats with Number=G, append _ref, _het, _hom. For things with Number=A, append _ref, _alt
    Specify which sample with its name or index in the VCF
    """
    v = pysam.VariantFile(fn)
    header = ["key", "id", "svtype", "svlen", "szbin", "qual", "filter", "is_pass"]
    fmts = []
    if with_fmt: # get all the format fields, and how to parse them from header, add them to the header
        for key in v.header.formats:
            num = v.header.formats[key].number
            if num in [1, '.', "A", 0]:
                header.append(key)
                fmts.append((key, lambda x: [x] if x is not None else [None]))
            elif num == 'R':
                header.append(key + '_ref')
                header.append(key + '_alt')
                fmts.append((key, lambda x: x if x is not None else [None, None]))
            elif num == 'G':
                header.append(key + '_ref')
                header.append(key + '_het')
                header.append(key + '_hom')
                fmts.append((key, lambda x: x if x is not None else [None, None, None]))
            else:
                logging.critical("Unknown Number (%s) for %s. Skipping.", num, key)

    infos = list(v.header.info.keys()) if with_info else []
    header.extend(infos)

    rows = []
    for entry in v:
        varsize = truvari.entry_size(entry)
        filt = list(entry.filter)
        cur_row = [f"{entry.chrom}:{entry.start}-{entry.stop}.{entry.alts[0]}",
                    entry.id,
                    truvari.entry_variant_type(entry),
                    varsize,
                    truvari.get_sizebin(varsize),
                    entry.qual,
                    filt,
                    filt == [] or filt[0] == 'PASS'
                  ]
        for i, op in fmts:
            if i in entry.samples[sample]:
                cur_row.extend(op(entry.samples[sample][i]))
            else:
                cur_row.extend(op(None))

        for i in infos:
            if i in entry.info:
                cur_row.append(entry.info[i])
            else:
                cur_row.append(None)
        rows.append(cur_row)
    ret = pd.DataFrame(rows, columns=header)
    ret["szbin"] = ret["szbin"].astype(SZBINTYPE)
    ret["svtype"] = ret["svtype"].astype(SVTYTYPE)
    return ret.set_index("key")

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="vcf2df", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf", metavar="VCF",
                        help="VCF to parse")
    parser.add_argument("output", metavar="JL",
                        help="Output joblib to save")
    parser.add_argument("-b", "--bench-dir", action="store_true",
                        help="Input is a truvari bench directory")
    parser.add_argument("-i", "--info", action="store_true",
                        help="Attempt to put the INFO fields into the dataframe")
    parser.add_argument("-f", "--format", action="store_true",
                        help="Attempt to put the FORMAT fileds into the dataframe")
    parser.add_argument("-S", "--skip-compression", action="store_true",
                        help="Skip the attempt to optimize the dataframe's size")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    return args

def vcf2df_main(args):
    """
    Main entry point for running DataFrame builder
    """
    args = parse_args(args)

    out = None
    if not args.bench_dir:
        out = vcf_to_df(args.vcf, args.info, args.format)
    else:
        vcfs = get_files_from_truvdir(args.vcf)
        all_dfs = []
        for key, val in vcfs.items():
            df = vcf_to_df(val[0], args.info, args.format)
            df["state"] = key
            all_dfs.append(df)
        out = pd.concat(all_dfs)

    # compression -- this is not super important for most VCFs
    # Especially since nulls are highly likely
    logging.info("Optimizing memory")
    pre_size = out.memory_usage().sum()
    for col in out.columns:
        try:
            if out[col].apply(float.is_integer).all():
                if len(out[out[col] < 0]) == 0:
                    out[col] = pd.to_numeric(out[col], downcast="unsigned")
                else:
                    out[col] = pd.to_numeric(out[col], downcast="signed")
            else:
                out[col] = pd.to_numeric(out[col], downcast="float")
        except TypeError as e:
            logging.debug("Unable to downsize %s (%s)", col, str(e))
    post_size = out.memory_usage().sum()
    logging.info("Optimized %.2fMB to %.2fMB", pre_size / 1e6, post_size / 1e6)
    joblib.dump(out, args.output)
    logging.info("Finished")


if __name__ == '__main__':
    vcf2df_main(sys.argv[1:])
