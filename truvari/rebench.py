"""
Automated Truvari bench result refinement
"""
import os
import sys
import json
import logging
import argparse
from argparse import Namespace
import itertools
from dataclasses import dataclass

import pysam
from intervaltree import IntervalTree

import truvari
from truvari.bench import StatsBox, compare_chunk, setup_outputs, output_writer, close_outputs

def read_json(fn):
    """
    Parse json and return dict
    """
    ret = None
    with open(fn, 'r') as fh:
        ret = json.load(fh)
    return ret

def intersect_beds(bed_a, bed_b):
    """
    Return bed_a regions that intersect bed_b regions
    """
    shared = {}
    count = 0
    for chrom in bed_a:
        s_inc = []
        for i in bed_a[chrom]:
            if bed_b[chrom].overlaps(i):
                count += 1
                s_inc.append(i)
        shared[chrom] = IntervalTree(s_inc)
    return shared

def update_fns(benchdir, fn_trees, summary):
    """
    For all FP calls in the fn_trees regions, update the FN counts
    summary is updated in-place
    """
    vcf = pysam.VariantFile(os.path.join(benchdir, "fn.vcf.gz"))
    for chrom in fn_trees:
        for intv in fn_trees[chrom]:
            for _ in vcf.fetch(chrom, intv.begin, intv.end):
                summary["FN"] -= 1
                summary["base cnt"] -= 1
@dataclass
class ReevalRegion: # pylint: disable=too-many-instance-attributes
    """
    Class for keeping track of a region for processing.
    """
    name: str
    benchdir: str
    params: dict
    reference: str
    chrom: str
    start: int
    end: int
    eval_method: str = "phab"
    needs_reeval: bool = True
    # These may need to become a list
    #base_calls: List = field(default_factory=lambda: [])
    base_calls: str = None
    base_count: int = 0
    #comp_calls: List = field(default_factory=lambda: [])
    comp_calls: str = None
    call_count: int = 0
    in_tp_base_count: int = 0
    in_tp_call_count: int = 0
    in_fp_count: int = 0
    in_fn_count: int = 0
    out_tp_base_count: int = 0
    out_tp_call_count: int = 0
    out_fp_count: int = 0
    out_fn_count: int = 0

    def set_out_to_in(self):
        """
        For regions with no changes, the output counts are identical to input counts
        """
        self.out_tp_base_count = self.in_tp_base_count
        self.out_tp_call_count = self.in_tp_call_count
        self.out_fp_count = self.in_fp_count
        self.out_fn_count = self.in_fn_count

    def update_summary(self, summary):
        """
        Given a StatsBox (dict), update the counts based on what we've observed

        This is currently done stupidly. I'm sure there's a cleaner logic for updating the counts,
        but I'll refactor it later
        """
        summary["TP-base"] += self.out_tp_base_count
        summary["TP-call"] += self.out_tp_call_count
        summary["FN"] += self.out_fn_count
        summary["FP"] += self.out_fp_count
        summary["base cnt"] += self.out_tp_base_count + self.out_fn_count
        summary["call cnt"] += self.out_tp_call_count + self.out_fp_count

    @staticmethod
    def get_header():
        """
        Get the header of the columns this makes when you output a string
        """
        return "\t".join(["chrom", "start", "end", "reevaled", "base_count", "call_count",
                          "in_tp_base_count", "in_tp_call_count", "in_fp_count", "in_fn_count",
                          "out_tp_base_count", "out_tp_call_count", "out_fp_count", "out_fn_count"])

    def __str__(self):
        return (f"{self.chrom}\t{self.start}\t{self.end}\t{self.needs_reeval}\t"
                f"{self.base_count}\t{self.call_count}\t{self.in_tp_base_count}\t"
                f"{self.in_tp_call_count}\t{self.in_fp_count}\t{self.in_fn_count}\t"
                f"{self.out_tp_base_count}\t{self.out_tp_call_count}\t{self.out_fp_count}\t{self.out_fn_count}")

def region_needs_reeval(region):
    """
    Determines if a ReevalRegion even needs to be re-evaluated

    Note: eventually I would like this method to also determine which EVALS method is best.
    When we implement that, we'll add a new member to the data class and... I might can
    design this in now...
    """
    region.needs_reeval = region.in_fn_count > 0 and region.in_fp_count > 0
    if not region.needs_reeval:
        region.set_out_to_in()

    return region

def pull_variants(in_vcf_fn, out_vcf_fn, chrom, start, end):
    """
    Move variants contained entirely inside region from in_vcf to out_vcf

    Returns how many variants were pulled
    """
    in_vcf = pysam.VariantFile(in_vcf_fn)
    out_vcf = pysam.VariantFile(out_vcf_fn, mode='w', header=in_vcf.header)
    cnt = 0
    for entry in in_vcf.fetch(chrom, start, end):
        ent_start, ent_end = truvari.entry_boundaries(entry)
        if start <= ent_start and ent_end <= end:
            cnt += 1
            out_vcf.write(entry)
    out_vcf.close()
    truvari.compress_index_vcf(out_vcf_fn)
    return cnt

def fetch_originals(region):
    """
    Grabs the variants that are needed for processing

    Populates the RevalRegion's `base_calls` and `comp_calls`

    ToDo: This should also check that the calls are inside the boundary
    """
    fetch_bench_vcfs(region, populate_calls=False)
    bout_name = truvari.make_temp_filename(suffix=".vcf")
    region.base_count = pull_variants(region.params["base"], bout_name, region.chrom, region.start, region.end)
    region.base_calls = bout_name + '.gz'

    cout_name = truvari.make_temp_filename(suffix=".vcf")
    region.call_count = pull_variants(region.params["comp"], cout_name, region.chrom, region.start, region.end)
    region.comp_calls = cout_name + '.gz'

    # if/when we do in-memory parsing of calls
    #if region.eval_method == 'phab':
    #else:
    # region.base_calls = [_ for _ in pysam.VariantFile(region.params["base"]).fetch(region.chrom, region.start, region.end)]
    # region.comp_calls = [_ for _ in pysam.VariantFile(region.params["comp"]).fetch(region.chrom, region.start, region.end)]
    return region

def fetch_bench_vcfs(region, populate_calls=True):
    """
    Grabs the variants that are needed for processing
    Populate the input state counts from truvari vcfs

    Populates the RevalRegion's `base_calls` and `comp_calls`

    ToDo: Clean this up.. it's gotta be simpler than this. I already tried to make pull_variants work for me, but it
    doesn't exactly work. So refactor that
    """
    tp_vcf = pysam.VariantFile(os.path.join(region.benchdir, "tp-base.vcf.gz"))
    bout_name, bout, cout_name, cout = None, None, None, None
    b_count, c_count = 0, 0
    if populate_calls:
        bout_name = truvari.make_temp_filename(suffix=".vcf")
        bout = pysam.VariantFile(bout_name, 'w', header=tp_vcf.header)

    for entry in tp_vcf.fetch(region.chrom, region.start, region.end):
        region.in_tp_base_count += 1
        b_count += 1
        if populate_calls:
            bout.write(entry)
            #region.base_calls.append(entry)

    fn_vcf = pysam.VariantFile(os.path.join(region.benchdir, "fn.vcf.gz"))
    for entry in fn_vcf.fetch(region.chrom, region.start, region.end):
        region.in_fn_count += 1
        b_count += 1
        if populate_calls:
            bout.write(entry)
            #region.base_calls.append(entry)

    tpc_vcf = pysam.VariantFile(os.path.join(region.benchdir, "tp-call.vcf.gz"))
    if populate_calls:
        cout_name = truvari.make_temp_filename(suffix=".vcf")
        cout = pysam.VariantFile(cout_name, 'w', header=tpc_vcf.header)

    for entry in tpc_vcf.fetch(region.chrom, region.start, region.end):
        region.in_tp_call_count += 1
        c_count += 1
        if populate_calls:
            cout.write(entry)
            #region.comp_calls.append(entry)

    fp_vcf = pysam.VariantFile(os.path.join(region.benchdir, "fp.vcf.gz"))
    for entry in fp_vcf.fetch(region.chrom, region.start, region.end):
        region.in_fp_count += 1
        c_count += 1
        if populate_calls:
            cout.write(entry)
            #region.comp_calls.append(entry)

    if populate_calls:
        bout.close()
        cout.close()
        truvari.compress_index_vcf(bout_name)
        truvari.compress_index_vcf(cout_name)
        region.base_calls = bout_name + '.gz'
        region.base_count = b_count
        region.comp_calls = cout_name + '.gz'
        region.call_count = c_count

    return region

def phab_eval(region):
    """
    Analyze a ReevalRegion with phab
    """
    if not region.needs_reeval:
        return region

    phab_dir = os.path.join(region.benchdir, "phab", region.name)
    os.makedirs(phab_dir)
    m_reg = (region.chrom, region.start, region.end)
    truvari.phab(region.base_calls, region.reference, phab_dir, m_reg,
                 comp_vcf=region.comp_calls, prefix_comp=True)

    # Setup/run truvari bench
    vcf_fn = os.path.join(phab_dir, "output.vcf.gz")
    m_args = Namespace(**region.params)
    m_args.no_ref = 'a'
    m_args.output = os.path.join(phab_dir, "bench")
    m_args.base = vcf_fn
    m_args.comp = vcf_fn
    outputs = run_bench(m_args)

    # Refine counts
    box = outputs["stats_box"]
    region.out_fn_count = box["FN"]
    region.out_fp_count = box["FP"]
    region.out_tp_base_count = box["TP-base"]
    region.out_tp_call_count = box["TP-call"]

    return region

def hap_eval(region):
    """
    Analyze a ReevalRegion with hap_eval
    """
    # Short circuit
    if not region.needs_reeval:
        return region
    raise NotImplementedError("In progress")

def run_bench(m_args):
    """
    Run truvari bench.

    Returns the bench outputs dict built by truvari.setup_outputs

    For now takes a single parameter `m_args` - a Namespace of all the command line arguments
    used by bench.

    This puts the burden on the user to
        1. build that namespace correctly (there's no checks on it)
        2. know how to use that namespace to get their pre-saved vcf(s) through
        3. read/process the output vcfs
        4. understand the `setup_outputs` structure even though that isn't an object

    Future versions I'll clean this up to not rely on files. Would be nice to have a way to just provide
    lists of base/comp calls and to return the e.g. output vcf entries with an in-memory object(s)

    This current version is just a quick convience thing

    Even with this quick thing which is almost essentially a command line wrapper, I could make it better
    with:
        make a bench params dataclass - helps document/standardize m_args
        make a `setup_outputs` dataclass - helps document/standardize outputs
    """
    matcher = truvari.Matcher(args=m_args)
    outputs = setup_outputs(m_args, do_logging=False)
    base = pysam.VariantFile(m_args.base)
    comp = pysam.VariantFile(m_args.comp)
    regions = truvari.RegionVCFIterator(base, comp, max_span=m_args.sizemax)
    base_i = regions.iterate(base)
    comp_i = regions.iterate(comp)
    chunks = truvari.chunker(matcher, ('base', base_i), ('comp', comp_i))
    for match in itertools.chain.from_iterable(map(compare_chunk, chunks)):
        output_writer(match, outputs, m_args.sizemin)
    box = outputs["stats_box"]
    with open(os.path.join(m_args.output, "summary.json"), 'w') as fout:
        box.calc_performance()
        fout.write(json.dumps(box, indent=4))
        logging.debug("%s Stats: %s", m_args.output, json.dumps(box, indent=4))

    close_outputs(outputs, True)
    return outputs


EVALS = {'phab':phab_eval, 'hap':hap_eval}

def rebench(benchdir, params, summary, reeval_trees, reference, use_original=False, eval_method='phab', workers=1):
    """
    Sets up inputs before calling specified EVAL method
    Will take the output from the EVAL method to then update summary

    :param `benchdir`: Truvari bench output directory
    :type `benchdir`: :class:`str`
    :param `params`: bench params.json
    :type `benchdir`: :class:`dict`
    :param `summary`: bench summary.json
    :type `bench`: :class:`dict`
    :param `reeval_trees`: dict of IntervalTrees of regions to be re-evaluated
    :type `reeval_trees`: `dict`
    :param `reference`: Reference to use for the eval method
    :type `reference`: str
    :param `use_original`: true if use original VCFs, else parse calls in benchdir VCFs
    :type `use_original`: `bool`
    :param `eval_method`:
    :param `workers`: number of workers to use
    :type `workers`: `int`
    """
    # Build a Region datastructure.
    # These are what we use to pass inputs/output
    data = []
    for chrom in reeval_trees:
        for intv in reeval_trees[chrom]:
            name = f"{chrom}:{intv.begin}-{intv.end}"
            data.append(ReevalRegion(name, benchdir, params, reference, chrom, intv.begin, intv.end, eval_method))

    # Build the pipeline
    pipeline = []

    # Get/count input variants
    if use_original:
        pipeline.append((fetch_originals))
    else:
        pipeline.append((fetch_bench_vcfs))

    # Is this even a region worth considering?
    pipeline.append((region_needs_reeval))

    # Call the eval_method
    pipeline.append((EVALS[eval_method]))

    # Collect results
    with open(os.path.join(benchdir, 'rebench.counts.txt'), 'w') as fout:
        fout.write(ReevalRegion.get_header() + '\n')
        for result in truvari.fchain(pipeline, data, workers=workers):
            result.update_summary(summary)
            fout.write(str(result) + '\n')

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="rebench", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("benchdir", metavar="DIR",
                        help="Truvari bench directory")
    parser.add_argument("-f", "--reference", type=str,
                        help="Indexed fasta used to call variants")
    parser.add_argument("-r", "--regions", default=None,
                        help="Regions to process")
    # commenting out until we actually have other eval methods
    #parser.add_argument("-e", "--eval", default='phab', choices=EVALS.keys(),
    #                    help="Evaluation procedure (%(default)s)")
    parser.add_argument("-u", "--use-original", action="store_true",
                        help="Use original input VCFs instead of filtered tp/fn/fps")
    parser.add_argument("-t", "--threads", default=1, type=int,
                        help="Number of threads to use (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug, show_version=True)
    return args

def rebench_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)
    param_path = os.path.join(args.benchdir, "params.json")
    if not os.path.exists(param_path):
        logging.error("Bench directory %s doesn't have params.json", param_path)
        sys.exit(1)
    params = read_json(param_path)

    if params["reference"] is None and args.reference is None:
        logging.error("Reference not in params.json or given as a parameter to rebench")
        sys.exit(1)
    elif args.reference is None:
        args.reference = params["reference"]

    if not os.path.exists(args.reference):
        logging.error("Reference %s does not exist", args.reference)
        sys.exit(1)

    # Setup prefix
    params["cSample"] = "p:" + params["cSample"]

    summary_path = os.path.join(args.benchdir, "summary.json")
    if not os.path.exists(summary_path):
        logging.error("Bench directory %s doesn't have summary.json", param_path)
        sys.exit(1)

    # Should check that bench dir has compressed/indexed vcfs for fetching

    #if args.eval == 'phab':
    #    # Might exist.. so another thing we gotta check
    phdir = os.path.join(args.benchdir, 'phab')
    if os.path.exists(phdir):
        logging.error("Directory %s exists. Cannot run phab", phdir)
        sys.exit(1)
    os.makedirs(phdir)

    if params["includebed"] is None and args.regions is None:
        logging.error("Bench output didn't use `--includebed` and `--regions` not provided")
        logging.error("Unable to run rebench")
        sys.exit(1)
    elif args.regions is None:
        reeval_trees, new_count = truvari.build_anno_tree(params["includebed"], idxfmt="")
        logging.info("rebench-ing %d regions", new_count)
    elif args.regions is not None and params["includebed"] is not None:
        a_trees, regi_count = truvari.build_anno_tree(args.regions, idxfmt="")
        b_trees, orig_count = truvari.build_anno_tree(params["includebed"], idxfmt="")
        reeval_trees, new_count  = intersect_beds(a_trees, b_trees)
        logging.info("%d --regions reduced to %d after intersecting with %d from --includebed",
                     regi_count, new_count, orig_count)
        # might need to merge overlaps.

    #update_fns(args.benchdir, fn_trees, summary) - no longer do this

    # Will eventually need to pass args for phab and|or hap-eval
    summary = StatsBox()
    rebench(args.benchdir, params, summary, reeval_trees,
            args.reference, args.use_original, 'phab', args.threads)
    summary.calc_performance()

    with open(os.path.join(args.benchdir, 'rebench.summary.json'), 'w') as fout:
        json.dump(summary, fout, indent=4)
    logging.info(json.dumps(summary, indent=4))
    logging.info("Finished rebench")
