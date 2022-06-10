"""
Over multiple vcfs, calculate their intersection/consistency.

Calls will match between VCFs if they have a matching key of:
    CHROM:POS ID REF ALT
"""
# pylint: disable=consider-using-f-string
import gzip
import argparse

from collections import defaultdict, Counter


def parse_vcf(fn):
    """
    Simple vcf reader
    """
    openfn = gzip.open if fn.endswith(".gz") else open
    with openfn(fn, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            # Only keep the first 5 fields, and use them as key
            yield "\t".join(line.split("\t")[:5])


def read_files(allVCFs):
    """
    Read all VCFs and mark all (union) calls for their presence

    For each call, we will use and integer and mark a call's presence
    at the index of the VCF file as 1.

    For example, if we have 3 VCFs, and the call is present in only the
    first VCF, then the integer will be `0 | (1 << 2)`, that is `0b100`,
    where the first bit marks that the call appears in the first VCF.

    Returns:
        all_presence: dict of all calls, with the integers as values.
            The integers are the bitwise OR of all VCFs that have the call.
        n_calls_per_vcf: list of the number of calls in each VCF
    """
    n_vcfs = len(allVCFs)
    # Initialize the integer to 0
    all_presence = defaultdict(lambda: 0)
    n_calls_per_vcf = []
    for i, vcf in enumerate(allVCFs):
        n_calls_per_vcf.append(0)
        for key in parse_vcf(vcf):
            n_calls_per_vcf[i] += 1
            all_presence[key] |= (1 << (n_vcfs - i - 1))

    # We don't care about the calls anyway for stats
    # Then this becomes a list of integers, which will save memory
    return all_presence.values(), n_calls_per_vcf


def get_shared_calls(all_presence, n):
    """Get n shared calls from the all_presence dictionary"""
    return sum(
        1 for presence in all_presence
        if bin(presence).count("1") == n
    )


def parse_args(args):
    """ parse args """
    parser = argparse.ArgumentParser(prog="consistency", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("allVCFs", metavar='VCFs', nargs='+',
                        help="VCFs to intersect")
    args = parser.parse_args(args)
    return args


def write_report(allVCFs, all_presence, n_calls_per_vcf):
    """
    Write the report
    """
    total_unique_calls = len(all_presence)
    n_vcfs = len(allVCFs)

    print("#\n# Total %d calls across %d VCFs\n#" %
          (total_unique_calls, n_vcfs))
    print("#File\tNumCalls")
    for i, fn in enumerate(allVCFs):
        print("%s\t%d" % (fn, n_calls_per_vcf[i]))

    print("#\n# Summary of consistency\n#")
    print("#VCFs\tCalls\tPct")

    for i in reversed(range(n_vcfs)):
        shared_calls = get_shared_calls(all_presence, i + 1)
        print("%d\t%d\t%.2f%%" % (
            i + 1,
            shared_calls,
            shared_calls / total_unique_calls * 100
        ))

    print("#\n# Breakdown of VCFs' consistency\n#")
    print("#Group\tTotal\tTotalPct\tPctOfFileCalls")

    all_overlap = Counter(all_presence)
    for group, ncalls in sorted(all_overlap.items(), key=lambda x: (-x[1], x[0])):
        # '0b100' -> '100'
        group = bin(group)[2:].rjust(n_vcfs, '0')
        c_text = " ".join(
            "%.2f%%" % (ncalls / n_calls_per_vcf[i] * 100)
            if group[i] == "1"
            else "0%"
            for i in range(n_vcfs)
        )
        print("%s\t%d\t%.2f%%\t%s" %
              (group, ncalls, ncalls / total_unique_calls * 100, c_text))


def consistency_main(args):
    """
    Run the program
    """
    args = parse_args(args)

    all_presence, n_calls_per_vcf = read_files(args.allVCFs)

    write_report(args.allVCFs, all_presence, n_calls_per_vcf)
