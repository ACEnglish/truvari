"""
Takes a vcf and creates a data frame. Can parse a bench output directory
"""
import os
import sys
import glob
import logging
import argparse
from enum import Enum

import pysam
import joblib
import pandas as pd
import truvari


class GT(Enum):
    """Genotypes

    - HET = <GT.HET: 1>
    - HOM = <GT.HOM: 2>
    - NON = <GT.NON: 3> - Non-Genotyped (e.g. `./.`)
    - REF = <GT.REF: 0>
    - UNK = <GT.UNK: 4> - Undetermined Genotype
    """
    NON = 3 
    REF = 0
    HET = 1
    HOM = 2
    UNK = 4


class SV(Enum):
    """SVtypes

    - DEL = <GT.HET: 0>
    - INS = <GT.HOM: 1>
    - DUP = <GT.HOM: 2>
    - INV = <GT.HOM: 3>
    - NON = <GT.HOM: 4> - Not an SV, SVTYPE
    - UNK = <GT.HOM: 5> - Unknown SVTYPE
    """
    DEL = 0
    INS = 1
    DUP = 2
    INV = 3
    NON = 4  # Not an SV, SVTYPE
    UNK = 5  # Unknown SVTYPE


SZBINS = ["[0,50)", "[50,100)", "[100,200)", "[200,300)", "[300,400)",
          "[400,600)", "[600,800)", "[800,1k)", "[1k,2.5k)",
          "[2.5k,5k)", ">=5k"]
SZBINMAX = [50, 100, 200, 300, 400, 600, 800, 1000, 2500, 5000, sys.maxsize]
QUALBINS = [f"[{x},{x+10})" for x in range(0, 100, 10)] + [">=100"]
SZBINTYPE = pd.CategoricalDtype(categories=SZBINS, ordered=True)
SVTYTYPE = pd.CategoricalDtype(categories=[_.name for _ in SV], ordered=True)


def get_svtype(svtype):
    """
    Turn to an SV

    :param `svtype`: SVTYPE string to turn into SV object
    :type `svtype`: string

    :return: A :class:`SV` of the SVTYPE
    :rtype: :class:`truvari.SV`
    """
    try:
        return SV.__members__[svtype]
    except AttributeError:
        pass
    return SV.UNK


def get_sizebin(sz):
    """
    Bin a given size into :data:`truvari.SZBINS`

    :param `sz`: SVLEN to bin into the SZBINS
    :type `sz`: integer

    :return: SZBIN
    :rtype: string
    """
    sz = abs(sz)
    for key, maxval in zip(SZBINS, SZBINMAX):
        if sz < maxval:
            return key
    return None


def get_gt(gt):
    """
    Turn a genotype tuple into a GT object

    :param `gt`: Genotype tuple
    :type `gt`: tuple (int, int)

    :return: A :class:`GT` of the genotype
    :rtype: :class:`truvari.GT`
    """
    if None in gt:
        return GT.NON
    if len(gt) != 2:
        return GT.UNK
    if gt == (0, 0):
        return GT.REF
    if gt[0] != gt[1]:
        return GT.HET
    if gt[0] == gt[1]:
        return GT.HOM
    return GT.UNK


def get_scalebin(x, rmin=0, rmax=100, tmin=0, tmax=100, step=10):
    """
    Scale variable x from rdomain to tdomain with step sizes
    return key, index

    :param `x`: Number to scale
    :type `x`: number
    :param `rmin`: The minimum of the range of your measurement
    :type `rmin`: number, optional
    :param `rmax`: The maximum of the range of your measurement
    :type `rmax`: number, optional
    :param `tmin`: The minimum of the range of your desired target scaling
    :type `tmin`: number, optional
    :param `tmax`: The maximum of the range of your measurement
    :type `tmax`: number, optional
    :param `step`: The step size of bins of target range
    :type `step`: number, optional

    :return: The bin-string of the scaled variable
    :rtype: string

    >>> import truvari
    >>> truvari.get_scalebin(4, 1, 5, 0, 20, 5)
    ('[15,20)', 3)
    """
    newx = (x - rmin) / (rmax - rmin) * (tmax - tmin) + tmin
    pos = 0
    for pos, i in enumerate(range(tmin, tmax, step)):
        if newx < i + step:
            return f"[{i},{i+step})", pos
    return f">={tmax}", pos + 1


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
            logging.error(
                "Expected 1 file for %s[.gz] found %d", fname, len(ret[key]))
            has_error = True
    #summary = glob.glob(os.path.join(directory, "summary.txt"))
    # ret["summary"]
    if has_error:
        raise FileNotFoundError(
            f"Couldn't parse truvari directory {directory}. See above errors")
    return ret


def vcf_to_df(fn, with_info=True, with_fmt=True, sample=0):
    """
    Parse a vcf file and turn it into a dataframe.
    Tries its best to pull info/format tags from the sample information
    For Formats with Number=G, append _ref, _het, _hom. For things with Number=A, append _ref, _alt
    Specify which sample with its name or index in the VCF

    :param `fn`: File name of VCF to open and turn into a DataFrame
    :type `fn`: string
    :param `with_info`:  Add the INFO fields from the VCF to the DataFrame columns
    :type `with_info`: boolean, optional
    :param `with_fmt`: Add the FORMAT fields from the VCF to the DataFrame columns
    :type `with_info`: boolean, optional
    :param `sample`: Sample from the VCF to parse. Only used when with_fmt==True
    :type `sample`: int/string, optional

    :return: Converted VCF
    :rtype: pandas.DataFrame
    """
    v = pysam.VariantFile(fn)
    header = ["key", "id", "svtype", "svlen",
              "szbin", "qual", "filter", "is_pass"]
    fmts = []
    if with_fmt:  # get all the format fields, and how to parse them from header, add them to the header
        for key in v.header.formats:
            num = v.header.formats[key].number
            if num in [1, '.', "A", 0]:
                header.append(key)
                fmts.append((key, lambda x: [x] if x is not None else [None]))
            elif num == 'R':
                header.append(key + '_ref')
                header.append(key + '_alt')
                fmts.append(
                    (key, lambda x: x if x is not None else [None, None]))
            elif num == 'G':
                header.append(key + '_ref')
                header.append(key + '_het')
                header.append(key + '_hom')
                fmts.append(
                    (key, lambda x: x if x is not None else [None, None, None]))
            else:
                logging.critical(
                    "Unknown Number (%s) for %s. Skipping.", num, key)

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
    parser.add_argument("-s", "--sample", default=0,
                        help="SAMPLE name to parse when building columns for --format")
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
        out = vcf_to_df(args.vcf, args.info, args.format, args.sample)
    else:
        vcfs = get_files_from_truvdir(args.vcf)
        all_dfs = []
        for key, val in vcfs.items():
            df = vcf_to_df(val[0], args.info, args.format, args.sample)
            df["state"] = key
            all_dfs.append(df)
        out = pd.concat(all_dfs)
        state_type = pd.CategoricalDtype(
            categories=['tpbase', 'fn', 'tp', 'fp'], ordered=True)
        out["state"] = out["state"].astype(state_type)

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
    logging.info("Finished vcf2df")
