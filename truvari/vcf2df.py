"""
Takes a vcf and creates a data frame. Can parse a bench output directory
"""
import os
import sys
import glob
import logging
import argparse
import itertools
from enum import Enum

import pysam
import joblib
import pandas as pd
import truvari


class GT(Enum):
    """Genotypes

    - REF = <GT.REF: 0>
    - HET = <GT.HET: 1>
    - HOM = <GT.HOM: 2>
    - NON = <GT.NON: 3> - Non-Genotyped (e.g. `./.`)
    - UNK = <GT.UNK: 4> - Undetermined Genotype
    """
    REF = 0
    HET = 1
    HOM = 2
    NON = 3
    UNK = 4


class SV(Enum):
    """SVtypes

    - SNP = <SV.SNP: 0>
    - DEL = <SV.DEL: 1>
    - INS = <SV.INS: 2>
    - DUP = <SV.DUP: 3>
    - INV = <SV.INV: 4>
    - NON = <SV.NON: 5> - Not a variant (monomorphic ref?)
    - UNK = <SV.UNK: 6> - Unknown type
    """
    SNP = 0
    DEL = 1
    INS = 2
    DUP = 3
    INV = 4
    NON = 5  # Not a variant (monomorphic ref?)
    UNK = 6  # Unknown type


SZBINS = ['SNP', '[1,5)', '[5,10)', '[10,15)', '[15,20)', '[20,30)', '[30,40)',
          '[40,50)', '[50,100)', '[100,200)', '[200,300)', '[300,400)',
          '[400,600)', '[600,800)', '[800,1k)', '[1k,2.5k)', '[2.5k,5k)', '>=5k']
SZBINMAX = [1, 5, 10, 15, 20, 30, 40, 50, 100, 200,
            300, 400, 600, 800, 1000, 2500, 5000, sys.maxsize]
QUALBINS = [f"[{x},{x+10})" for x in range(0, 100, 10)] + [">=100"]
SZBINTYPE = pd.CategoricalDtype(categories=SZBINS, ordered=True)
SVTYTYPE = pd.CategoricalDtype(categories=[_.name for _ in SV], ordered=True)


def get_svtype(svtype):
    """
    Turn a SVTYPE string to a :class:`truvari.SV` object

    :param `svtype`:
    :type `svtype`: string

    :return: A :class:`SV` of the SVTYPE
    :rtype: :class:`truvari.SV`

    Example
        >>> import truvari
        >>> truvari.get_svtype("INS")
        <SV.INS: 2>
        >>> truvari.get_svtype("foo")
        <SV.UNK: 6>
    """
    try:
        return SV.__members__[svtype]
    except KeyError:
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

    Example
        >>> import truvari
        >>> truvari.get_gt((0, 0))
        <GT.REF: 0>
        >>> truvari.get_gt((0, 1, 1))
        <GT.UNK: 4>
    """
    if None in gt:  # Any . makes it NON
        return GT.NON
    if len(gt) != 2:
        return GT.UNK
    if gt == (0, 0):
        return GT.REF
    if gt[0] != gt[1]:
        return GT.HET
    # Must be
    return GT.HOM


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

    Example
        >>> import truvari
        >>> truvari.get_scalebin(4, 1, 5, 0, 20, 5)
        ('[15,20)', 3)
        >>> truvari.get_scalebin(6, 1, 5, 0, 20, 5)
        ('>=20', 4)
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
    pats = [("tpbase", "tp-base.vcf.gz"), ("tp", "tp-comp.vcf.gz"),
            ("fp", "fp.vcf.gz"), ("fn", "fn.vcf.gz")]
    has_error = False
    for key, name in pats:
        fname = os.path.join(directory, name)
        ret[key] = glob.glob(fname)
        if len(ret[key]) != 1:
            logging.error(
                "Expected 1 file for %s found %d", fname, len(ret[key]))
            has_error = True
    if has_error:
        raise FileNotFoundError((f"Couldn't parse truvari directory {directory}. "
                                 "See above errors"))
    return ret


def tags_to_ops(items):
    """
    Given a list of  items from a :class:`pysam.VariantFile.header.[info|format].items`
    build column names and operators for parsing

    :param `items`: SVTYPE string to turn into SV object
    :type `items`: list of tuples, (str, :class:`pysam.libcbcf.VariantMetadata`)

    :return: list of column names, list of (key, operator)
    :rtype: tuple
    """
    def prod(base, suffixes):
        return [f"{base}_{s}" for s in suffixes]

    def pres_check(dat, key, multi=False):
        if key not in dat:
            return False
        if dat[key] is None:
            return False
        if multi and dat[key].count(None) == len(dat[key]):
            return False
        return True

    columns = []
    ops = []
    for key, dat in items:
        num = dat.number
        if num in [1, 'A', '.']:  # 'A' should just pull first item...
            columns.append(key)
            ops.append(
                (key, lambda dat, k: [dat[k]] if pres_check(dat, k) else [None]))
        elif num == 0:
            columns.append(key)
            ops.append((key, lambda dat, k: [k in dat]))
        elif num == 'R':
            columns.extend(prod(key, ['ref', 'alt']))
            ops.append(
                (key, lambda dat, k: dat[k] if pres_check(dat, k, True) else [None, None]))
        elif num == 'G':
            columns.extend(prod(key, ['ref', 'het', 'hom']))
            ops.append(
                (key, lambda dat, k: dat[k] if pres_check(dat, k, True) else [None, None, None]))
        elif isinstance(num, int):
            columns.append(key)
            ops.append(
                (key, lambda dat, k: [dat[k]] if pres_check(dat, k) else [()]))
        else:
            logging.critical("Unknown Number (%s) for %s. Skipping.", num, key)
    return columns, ops


def vcf_to_df(fn, with_info=True, with_format=True, sample=None, no_prefix=False, alleles=False):
    """
    Parse a vcf file and turn it into a dataframe.
    Tries its best to pull info/format tags from the sample information.
    For Formats with Number=G, append _ref, _het, _hom. For things with Number=R, append _ref, _alt.
    Specify which sample with its name or index in the VCF.

    :param `fn`: File name of VCF to open and turn into a DataFrame
    :type `fn`: string
    :param `with_info`:  Add the INFO fields from the VCF to the DataFrame columns
    :type `with_info`: boolean, optional
    :param `with_format`: Add the FORMAT fields from the VCF to the DataFrame columns
    :type `with_info`: boolean, optional
    :param `sample`: Sample from the VCF to parse. Only used when with_format==True
    :type `sample`: int/string, optional
    :param `no_prefix`: Don't prefix FORMAT fields with sample name
    :type `no_prefix`: boolean, optional
    :param `alleles`: Add column with allele sequences
    :type `alleles`: boolean, optional

    :return: Converted VCF
    :rtype: pandas.DataFrame

    Example
        >>> import truvari
        >>> df = truvari.vcf_to_df("repo_utils/test_files/variants/input2.vcf.gz", True, True)
        >>> df.columns
        Index(['chrom', 'start', 'end', 'id', 'svtype', 'svlen', 'szbin', 'qual',
               'filter', 'is_pass', 'QNAME', 'QSTART', 'QSTRAND', 'SVTYPE', 'SVLEN',
               'NA12878_GT', 'NA12878_PL_ref', 'NA12878_PL_het', 'NA12878_PL_hom',
               'NA12878_AD_ref', 'NA12878_AD_alt'],
              dtype='object')
    """
    v = pysam.VariantFile(fn)
    if with_format and not sample:
        sample = list(v.header.samples)
    if sample and len(sample) > 1 and no_prefix:
        raise TypeError("Multiple samples being pulled, must use prefix")

    header = ["hash", "chrom", "start", "end", "id", "svtype", "svlen",
              "szbin", "qual", "filter", "is_pass"]

    info_ops = []
    if with_info:
        info_header, info_ops = tags_to_ops(v.header.info.items())
        logging.debug(info_header)
        header.extend(info_header)

    fmt_ops = []
    if with_format:  # get all the format fields, and how to parse them from header, add them to the header
        fmt_header, fmt_ops = tags_to_ops(v.header.formats.items())
        logging.debug(fmt_header)
        if isinstance(sample, list) and not no_prefix:
            header.extend(
                [f'{s}_{f}' for s, f in itertools.product(sample, fmt_header)])
        else:
            header.extend(fmt_header)

        if sample is None:
            sample = v.header.samples[0]

        if not isinstance(sample, list):
            sample = [sample]
    else:
        sample = []
    if alleles:
        header.append("alleles")

    def _transform():
        """
        Yields the rows
        """
        for entry in v:
            varsize = truvari.entry_size(entry)
            cur_row = [truvari.entry_to_hash(entry),
                       entry.chrom,
                       entry.start,
                       entry.stop,
                       entry.id,
                       truvari.entry_variant_type(entry).name,
                       varsize,
                       truvari.get_sizebin(varsize),
                       entry.qual,
                       list(entry.filter),
                       not truvari.entry_is_filtered(entry)
                       ]

            for i, op in info_ops:
                cur_row.extend(op(entry.info, i))

            for samp in sample:
                for i, op in fmt_ops:
                    cur_row.extend(op(entry.samples[samp], i))

            if alleles:
                cur_row.append(entry.alleles)

            yield cur_row

    ret = pd.DataFrame(_transform(), columns=header)
    ret["szbin"] = ret["szbin"].astype(SZBINTYPE)
    ret["svtype"] = ret["svtype"].astype(SVTYTYPE)
    return ret.set_index("hash")


def optimize_df_memory(df):
    """
    Optimize DataFrame memory by trying to downcast every column in-place.
    Returns the pre/post size

    :param `df`: Dataframe to optimize
    :type `df`: :class:`pandas.DataFrame`

    :return: tuple (pre_size, post_size)
    :rtype: int, sizes in bytes
    """
    pre_size = df.memory_usage().sum()
    for col in df.columns:
        try:
            if df[col].apply(float.is_integer).all():
                if len(df[df[col] < 0]) == 0:
                    df[col] = pd.to_numeric(df[col], downcast="unsigned")
                else:
                    df[col] = pd.to_numeric(df[col], downcast="signed")
            else:
                df[col] = pd.to_numeric(df[col], downcast="float")
        except TypeError as e:
            logging.debug("Unable to downsize %s (%s)", col, str(e))
    post_size = df.memory_usage().sum()
    return pre_size, post_size


def bench_dir_to_df(dir_name, *args, **kwargs):
    """
    Parse a bench directory - same interface as vcf_to_df, except for first argument
    """
    vcfs = get_files_from_truvdir(dir_name)
    all_dfs = []
    for key, val in vcfs.items():
        df = vcf_to_df(val[0], *args, **kwargs)
        df["state"] = key
        all_dfs.append(df)
    out = pd.concat(all_dfs)
    state_type = pd.CategoricalDtype(
        categories=['tpbase', 'fn', 'tp', 'fp'], ordered=True)
    out["state"] = out["state"].astype(state_type)
    return out


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
    parser.add_argument("-s", "--sample", default=None,
                        help="SAMPLE name to parse when building columns for --format")
    parser.add_argument("-n", "--no-prefix", action="store_true",
                        help="Don't prepend sample name to format columns")
    parser.add_argument("-S", "--skip-compression", action="store_true",
                        help="Skip the attempt to optimize the dataframe's size")
    parser.add_argument("-c", "--compress", metavar="LVL", type=int, default=3,
                        help="Compression level for joblib 0-9 (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    if args.sample:
        args.sample = args.sample.split(',')

    truvari.setup_logging(args.debug, show_version=True)
    if args.compress not in range(10):
        logging.warning('--compress must be between 0-9. Setting to 3')
        args.compress = 3
    return args


def vcf2df_main(args):
    """
    Main entry point for running DataFrame builder
    """
    args = parse_args(args)

    out = None
    if args.bench_dir:
        out = bench_dir_to_df(args.vcf, args.info,
                              args.format, args.sample, args.no_prefix)
    else:
        out = vcf_to_df(args.vcf, args.info, args.format,
                        args.sample, args.no_prefix)

    # compression -- this is not super important for most VCFs
    logging.info("Optimizing DataFrame memory")
    pre_size, post_size = optimize_df_memory(out)
    logging.info("Optimized %.2fMB to %.2fMB", pre_size / 1e6, post_size / 1e6)

    joblib.dump(out, args.output, compress=args.compress)
    logging.info("Finished vcf2df")
