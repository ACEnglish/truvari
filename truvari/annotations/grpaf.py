"""
Add allele frequency annotations for subsets of samples
"""
import sys
import logging
import argparse

import pysam
import pandas as pd
import truvari


def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="grpaf", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="VCF to annotate ")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output filename (stdout)")
    parser.add_argument("-l", "--labels", type=str, required=True,
                        help="Tab-delimited file of sample and group")
    parser.add_argument("-t", "--tags", type=str, default='all',
                        help=("Comma-separated list of tags to add "
                        "from AF,MAF,ExcHet,HWE,MAC,AC,AN (%(default)s)"))
    parser.add_argument("--strict", action="store_true",
                        help="Exit if sample listed in labels is not present in VCF (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    return args


def edit_header(header, tags, groups):
    """
    Add the header information
    """
    tmpl = '##INFO=<ID={mid},Type={mty},Number={mnum},Description="{desc}">'
    tag_meta = {}
    tag_meta["AF"] = ("Float", '1', "Allele Frequency on {count} {grp} samples")
    tag_meta["AN"] = ("Integer", '1', "Total number of alleles in called genotypes on {count} {grp} samples")
    tag_meta["MAF"] = ("Float", '1', "Minor Allele Frequency on {count} {grp} samples")
    tag_meta["AC"] = ("Integer", 'A', "Allele Count on {count} {grp} samples")
    tag_meta["MAC"] = ("Integer", 'A', "Minor Allele Count on {count} {grp} samples")
    tag_meta["HWE"] = ("Float", '1', "HWE test (PMID:15789306) on {count} {grp} samples; 1=good, 0=bad")
    tag_meta["ExcHet"] = ("Float", '1', "Test excess heterozygosity on {count} {grp} samples; 1=good, 0=bad")

    cnts = groups["group"].value_counts()
    for g_name, g_count in zip(cnts.index, cnts):
        for t in tags:
            m_id = f"{t}_{g_name}"
            m_ty, m_num, m_desc = tag_meta[t]
            m_desc = m_desc.format(grp=g_name, count=g_count)
            m_line = tmpl.format(mid=m_id, mty=m_ty, mnum=m_num, desc=m_desc)
            header.add_line(m_line)

    return header

def grpaf_main(cmd_args):
    """
    Main
    """
    args = parse_args(cmd_args)
    # validate tags
    all_tags = ["AF", "MAF", "ExcHet", "HWE", "MAC", "AC", "AN"]
    if args.tags == 'all':
        args.tags = all_tags
    else:
        args.tags = args.tags.split(',')
        unk_tags = set(args.tags) - set(all_tags)
        if unk_tags:
            logging.error("Unknown --tags %s", " ".join(unk_tags))
            logging.error("Use 'all' or %s", " ".join(all_tags))
            sys.exit(1)

    in_vcf = pysam.VariantFile(args.input)

    # setup masks
    vcf_samples = pd.Series(list(in_vcf.header.samples))
    groups = pd.read_csv(args.labels, sep='\t', names=["sample", "group"])
    presence = groups["sample"].isin(vcf_samples)
    if args.strict and not presence.all():
        logging.error("Samples in --labels not in VCF header")
        logging.error(f"\n{groups[~presence].to_string()}")
        sys.exit(1)
    group_lookup = {}
    for grp, vals in groups.groupby("group"):
        group_lookup[grp] = vcf_samples.isin(vals["sample"])

    # edit header
    n_header = edit_header(in_vcf.header.copy(), args.tags, groups)
    out = pysam.VariantFile(args.output, 'w', header=n_header)

    # process
    for entry in in_vcf:
        entry.translate(n_header)
        gt_series = pd.Series([_["GT"] for _ in entry.samples.values()])
        for grp, mask in group_lookup.items():
            af = truvari.calc_af(gt_series[mask])
            af["AC"] = af["AC"][1]
            for anno in args.tags:
                entry.info[f"{anno}_{grp}"] = af[anno]
        out.write(entry)

    logging.info("Finished grpaf")
