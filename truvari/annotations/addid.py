"""
Sets an ID to variants
"""
import logging
import argparse

import truvari


def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="addid", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", nargs="?", type=str, default="/dev/stdin",
                        help="VCF to annotate (stdin)")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output filename (stdout)")
    truvari.setup_logging(show_version=True)
    return parser.parse_args(args)


def get_idx(line):
    """
    Get string's tab positions
    """
    ret = [0]
    for _ in range(8):
        ret.append(line.index('\t', ret[-1] + 1))
    return ret


def addid_main(cmdargs):
    """
    Main method
    """
    args = parse_args(cmdargs)
    fh = truvari.opt_gz_open(args.input)

    cnt = 0
    with open(args.output, 'w') as fout:
        for line in fh:
            if line[0] == "#":
                fout.write(line)
                continue
            tabs = get_idx(line)
            new_id = line[:tabs[1]].replace('chr', '') + hex(cnt)
            fout.write(line[:tabs[2] + 1])
            fout.write(new_id)
            fout.write(line[tabs[3]:])
            cnt += 1
    logging.info("Finished addid")
