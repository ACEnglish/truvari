"""
Divide a VCF into independent parts
"""
import os
import sys
import logging
import argparse
import concurrent.futures

import pysam
import pandas as pd

import truvari

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="divide", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", metavar="VCF",
                        help="VCF to split")
    parser.add_argument("output", metavar="DIR",
                        help="Output directory to save parts")
    parser.add_argument("-b", "--buffer", type=truvari.restricted_int, default=1000,
                        help="Buffer to make mini-clusters (%(default)s)")
    parser.add_argument("-m", "--min", type=truvari.restricted_int, default=100,
                        help="Minimum number of entries per-vcf (%(default)s)")
    parser.add_argument("--no-compress", action="store_false",
                        help="Don't attempt to compress/index sub-VCFs")
    parser.add_argument("-T", "--threads", type=int, default=1,
                        help="Number of threads for compressing/indexing sub-VCFs (%(default)s)")
    args = parser.parse_args(args)
    truvari.setup_logging(False)
    return args

def flush_stack(in_vcf, stack, out_name, compress=True):
    """
    write the stack to out_name
    """
    logging.debug(f"{out_name} has {len(stack)} entries")
    cur_out = pysam.VariantFile(out_name, mode='w', header=in_vcf.header)
    for i in stack:
        cur_out.write(i)
    cur_out.close()
    if compress:
        truvari.compress_index(out_name)

def divide_main(args):
    """
    Main
    """
    args = parse_args(args)

    if os.path.exists(args.output):
        logging.error("Output directory exists")
        sys.exit(1)
    os.mkdir(args.output)

    in_vcf = pysam.VariantFile(args.input)

    oname = list(os.path.splitext(os.path.basename(args.input)))
    if oname[-1] == ".gz":
        oname.pop(-1)
        oname = list(os.path.splitext(".".join(oname)))

    if oname[-1] != ".vcf":
        logging.error("Expected VCF input. Pysam may raise error")
    oname.pop(-1)
    oname = ".".join(oname)

    out_name_template = os.path.join(args.output, oname) + "_{part_num:08d}.vcf"
    cluster_counts = []
    stack = [next(in_vcf)]

    max_end = stack[0].stop

    with concurrent.futures.ThreadPoolExecutor(max_workers = args.threads) as executor:
        futures = {}
        for entry in in_vcf:
            if entry.chrom != stack[0].chrom and len(stack) >= args.min:
                m_name = out_name_template.format(part_num=len(cluster_counts))
                futures[executor.submit(flush_stack, in_vcf, stack, m_name, args.no_compress)] = m_name
                #flush_stack(in_vcf, stack, m_name, args.no_compress)
                cluster_counts.append(len(stack))
                max_end = entry.stop
                stack = [entry]
            elif entry.start >= max_end + args.buffer and len(stack) >= args.min:
                m_name = out_name_template.format(part_num=len(cluster_counts))
                futures[executor.submit(flush_stack, in_vcf, stack, m_name, args.no_compress)] = m_name
                #flush_stack(in_vcf, stack, m_name, args.no_compress)
                cluster_counts.append(len(stack))
                stack = [entry]
                max_end = entry.stop
            else:
                stack.append(entry)
                max_end = max(max_end, entry.stop)

        if stack:
            m_name = out_name_template.format(part_num=len(cluster_counts))
            futures[executor.submit(flush_stack, in_vcf, stack, m_name, args.no_compress)] = m_name
            #flush_stack(in_vcf, stack, m_name, args.no_compress)
            cluster_counts.append(len(stack))

        for future in concurrent.futures.as_completed(futures):
            m_result = futures[future]
            try:
                future.result()
            except Exception as exc: # pylint: disable=broad-except
                logging.error('%s generated an exception: %s', m_result, exc)

    cluster_counts = pd.Series(cluster_counts)

    logging.info(cluster_counts.describe())
