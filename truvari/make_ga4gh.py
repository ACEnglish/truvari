"""
Consolidate truvari by truth/query and annotate with GA4GH intermediates tags
"""
import os
import sys
import json
import logging
import argparse
from collections import defaultdict

import pysam
import pandas as pd
from intervaltree import IntervalTree

import truvari


def parse_args(args):
    """
    Argument parser
    """
    parser = argparse.ArgumentParser(prog="ga4gh", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", required=True,
                        help="Truvari result directory")
    parser.add_argument("-o", "--output", required=True,
                        help="Output prefix")
    parser.add_argument("-r", "--no-refine", action="store_true",
                        help="Don't pull from refined results")
    parser.add_argument("-w", "--write-phab", action="store_true",
                        help="Count/Write the phab variant representations")
    parser.add_argument("-b", "--buffer", type=truvari.restricted_int, default=0,
                        help="Amount of buffer used around refined regions (%(default)s)")
    args = parser.parse_args(args)
    return args



def build_tree(regions, buffer=0):
    """
    Build tree from regions
    """
    tree = defaultdict(IntervalTree)
    for _, i in regions.iterrows():
        tree[i['chrom']].addi(i['start'] - buffer, i['end'] + buffer + 1)
    for i in tree:
        tree[i].merge_overlaps()
    return tree


def get_truvari_filenames(in_dir):
    """
    Join in_dir to relevant filenames
    """
    return {'tp-base': os.path.join(in_dir, 'tp-base.vcf.gz'),
            'tp-comp': os.path.join(in_dir, 'tp-comp.vcf.gz'),
            'fn': os.path.join(in_dir, 'fn.vcf.gz'),
            'fp': os.path.join(in_dir, 'fp.vcf.gz'),
            'params': os.path.join(in_dir, 'params.json')}



def check_bench_dir(dirname):
    """
    Make sure all the files are there
    """
    check_fail = False
    b_files = get_truvari_filenames(dirname)
    for k, v in b_files.items():
        if not os.path.exists(v):
            logging.error("%s bench file doesn't exist", k)
            check_fail = True
    return check_fail


def check_args(args):
    """
    Ensure everything we're looking for is available
    return True if check failed
    """
    check_fail = False
    if not os.path.isdir(args.input):
        logging.error("input is not a directory")
        check_fail = True
    else:
        check_fail |= check_bench_dir(args.input)

    if os.path.exists(args.output + '.base.vcf.gz'):
        logging.error("%s.base.vcf.gz already exists", args.output)
        check_fail = True
    if os.path.exists(args.output + '.comp.vcf.gz'):
        logging.error("%s.comp.vcf.gz already exists", args.output)
        check_fail = True

    if args.no_refine:
        return check_fail

    pdir = os.path.join(args.input, 'phab_bench')
    if not os.path.exists(pdir):
        logging.error("phab_bench dir doesn't exist")
        check_fail = True
    else:
        check_fail |= check_bench_dir(pdir)

    return check_fail


def pull_new_header(vcf_fn, sample):
    """
    Make a new header from a template
    """
    vcf = pysam.VariantFile(vcf_fn)
    header = pysam.VariantHeader()
    for i in vcf.header.records:
        out = str(i)
        if not out.startswith("##fileformat"):
            header.add_line(str(i))
    header.add_sample(sample)
    header.add_line(('##INFO=<ID=IsRefined,Number=0,Type=Flag,'
                     'Description="True if variant was pulled from refinement">'))
    header.add_line(('##FORMAT=<ID=BD,Number=1,Type=String,'
                     'Description="Decision for call (TP/FP/FN/N)">'))
    header.add_line(('##FORMAT=<ID=BK,Number=1,Type=String,'
                     'Description="Sub-type for decision (match/mismatch type) lm variants had any counterpart to which it could compare">'))
    return header


class Ga4ghOutput():
    """
    Helper class for consolidating benchmark vcfs into a ga4gh output format
    """
    def __init__(self, in_base, in_comp, out_prefix, bSample, cSample, buffer=0):
        # You gotta do the header copy stuff..
        # And also add the new FORMAT tags
        self.bSample = bSample
        b_header = pull_new_header(in_base, bSample)
        self.base = truvari.VariantFile(
            out_prefix + '.base.vcf', 'w', header=b_header)

        self.cSample = cSample
        c_header = pull_new_header(in_comp, cSample)
        self.comp = truvari.VariantFile(
            out_prefix + '.comp.vcf', 'w', header=c_header)
        self.stats = {'TP-base': 0, 'TP-comp': 0, 'FP': 0, 'FN': 0, 'precision': None,
                      'recall': None, 'f1': None, 'base cnt': 0, 'comp cnt': 0}
        self.buffer = buffer

    def pull_annotate(self, vcf, out_vcf, f_sample, t_sample, bdkey, is_refined=False):
        """
        vcf to pull from
        output vcf to write to
        f_sample name to pull from
        t_sample name to pull into
        state annotation
        """
        cnt = 0
        for entry in vcf:
            n_entry = entry.move_record(out_vcf, [f_sample], [t_sample])
            n_entry.info["IsRefined"] = is_refined
            n_entry.samples[t_sample]["BD"] = bdkey[:2]
            if bdkey.startswith('T'):
                if "GTMatch" in entry.info and entry.info["GTMatch"] == 0:
                    n_entry.samples[0]["BK"] = 'gm'
                else:
                    n_entry.samples[0]["BK"] = 'am'
            else:
                if "TruScore" in entry.info:
                    n_entry.samples[0]["BK"] = 'lm'
                else:
                    n_entry.samples[0]["BK"] = '.'
            cnt += 1
            out_vcf.write(n_entry)
        return cnt

    def pull_from_dir(self, pull_dir, regions=None, within=True, force_tp=False, from_b=0, from_c=0, is_refined=False):
        """
        Pull tp/fp/fn files and annotate BK/BD tags and update stats
        force_state None to keep the state, T or F to force the state prefix
        """
        vcfs = get_truvari_filenames(pull_dir)
        if regions is None:
            def puller(x):
                return x
        else:
            refine_tree = build_tree(regions, self.buffer)
            def puller(x):
                return x.fetch_regions(refine_tree, inside=within)

        # Base
        vcf = truvari.VariantFile(vcfs['tp-base'])
        self.stats['TP-base'] += self.pull_annotate(puller(vcf),
                                                    self.base,
                                                    from_b,
                                                    self.bSample,
                                                    'TP',
                                                    is_refined)
        vcf = truvari.VariantFile(vcfs['fn'])
        state = "TP-base" if force_tp else 'FN'
        self.stats[state] += self.pull_annotate(puller(vcf),
                                                self.base,
                                                from_b,
                                                self.bSample,
                                                state,
                                                is_refined)

        # Comp
        vcf = truvari.VariantFile(vcfs['tp-comp'])
        self.stats['TP-comp'] += self.pull_annotate(puller(vcf),
                                                    self.comp,
                                                    from_c,
                                                    self.cSample,
                                                    'TP',
                                                    is_refined)
        vcf = truvari.VariantFile(vcfs['fp'])
        state = "TP-comp" if force_tp else 'FP'
        self.stats[state] += self.pull_annotate(puller(vcf),
                                                self.comp,
                                                from_c,
                                                self.cSample,
                                                state,
                                                is_refined)

    def calc_stats(self):
        """
        Calculate stats of variants we just wrote
        """
        p, r, f = truvari.performance_metrics(self.stats['TP-base'],
                                              self.stats['TP-comp'],
                                              self.stats['FN'],
                                              self.stats['FP'])
        self.stats['base cnt'] = self.stats['TP-base'] + self.stats['FN']
        self.stats['comp cnt'] = self.stats['TP-comp'] + self.stats['FP']
        self.stats['precision'] = p
        self.stats['recall'] = r
        self.stats['f1'] = f

    def close(self):
        """
        Close VCFs
        """
        self.base.close()
        self.comp.close()


def make_ga4gh(input_dir, out_prefix, pull_refine=False, write_phab=False, subset=False):
    """
    Grab all the files from a bench directory
    refined_regions are the set of regions which underwent refinement, if None, tries to pull candidate.refine.bed
    write_phab will pull the variants from the result/phab_bench instead of update the annotation based
    on the refined_regions

    If subset, counts only from the refine.regions.txt are used. Only possible when pull_refine == True
    """
    bench_vcfs = get_truvari_filenames(input_dir)

    with open(os.path.join(input_dir, "params.json")) as fh:
        params = json.load(fh)

    output = Ga4ghOutput(bench_vcfs['tp-base'], bench_vcfs['tp-comp'],
                         out_prefix, params['bSample'], params['cSample'])

    if not pull_refine:
        output.pull_from_dir(input_dir)
    else:
        regions = pd.read_csv(os.path.join(
            input_dir, "refine.regions.txt"), sep='\t')
        regions = regions[regions['refined']] # Only interested in this subset
        if write_phab:
            phab_dir = os.path.join(input_dir, 'phab_bench')
            output.pull_from_dir(phab_dir, regions, within=True,
                                 from_b=params['bSample'],
                                 from_c=params['cSample'],
                                 is_refined=True)
        else:
            mask = regions['state'] == 'TP'
            # Reannotate the state as TP
            output.pull_from_dir(input_dir, regions[mask], within=True, force_tp=True,
                                 from_b=params['bSample'],
                                 from_c=params['cSample'],
                                 is_refined=True)
            # Keep the original states
            output.pull_from_dir(input_dir, regions[~mask], within=True,
                                 from_b=params['bSample'],
                                 from_c=params['cSample'],
                                 is_refined=True)
        # Already handled everything touched by refine
        # Now grab the rest if we're not subsetting
        if not subset:
            output.pull_from_dir(input_dir, regions, within=False,
                                 from_b=params['bSample'],
                                 from_c=params['cSample'])

    output.close()
    truvari.compress_index_vcf(output.base.filename.decode())
    truvari.compress_index_vcf(output.comp.filename.decode())

    output.calc_stats()
    return output


def make_ga4gh_main(args):
    """
    Main entrypoint
    """
    args = parse_args(args)
    truvari.setup_logging()

    if check_args(args):
        logging.error("Unable to do conversion")
        sys.exit(1)

    logging.info("Consolidating VCFs")
    output = make_ga4gh(args.input,
                        args.output,
                        pull_refine=not args.no_refine,
                        write_phab=args.write_phab)
    logging.info("Stats: %s", json.dumps(output.stats, indent=4))
    with open(args.output + '.summary.json' ,'w') as fout:
        json.dump(output.stats, fout, indent=4)
    logging.info("Finished")
