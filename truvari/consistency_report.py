"""
Over multiple vcfs, calculate their intersection/consistency.

Calls will match between VCFs if they have a matching key of:
    CHROM:POS ID REF ALT
"""
# pylint: disable=consider-using-f-string
import io
import gzip
import argparse
import json
from collections import defaultdict, Counter


def parse_vcf(fn):
    """
    Simple vcf reader
    """
    if fn.endswith(".gz"):
        fh = io.TextIOWrapper(gzip.open(fn))
    else:
        fh = open(fn, 'r')  # pylint: disable=consider-using-with
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
    parser.add_argument("-j", "--json", action='store_true',
                        help="Output report in json format")
    parser.add_argument("allVCFs", metavar='VCFs', nargs='+',
                        help="VCFs to intersect")
    args = parser.parse_args(args)
    return args


def make_report(allVCFs, all_presence, n_calls_per_vcf):
    """
    Create the report
    """
    output = {}
    output['vcfs'] = allVCFs
    total_unique_calls = len(all_presence)
    output['total_calls'] = total_unique_calls
    n_vcfs = len(allVCFs)
    output['num_vcfs'] = n_vcfs
    output['vcf_counts'] = {}
    for i, fn in enumerate(allVCFs):
        output['vcf_counts'][fn] = n_calls_per_vcf[i]
    output['shared'] = []
    for i in reversed(range(n_vcfs)):
        shared_calls = get_shared_calls(all_presence, i + 1)
        call_pct = shared_calls / total_unique_calls
        output['shared'].append({'vcf_count': i + 1, 'num_calls': shared_calls, 'call_pct': call_pct})
    output['detailed'] = []
    all_overlap = Counter(all_presence)
    for group, ncalls in sorted(all_overlap.items(), key=lambda x: (-x[1], x[0])):
        group = bin(group)[2:].rjust(n_vcfs, '0')
        tot_pct = ncalls / total_unique_calls
        c = [(ncalls / n_calls_per_vcf[i])
             if group[i] == "1"
             else 0
             for i in range(n_vcfs)
        ]
        row = {'group': group,
               'total': ncalls,
               'total_pct': tot_pct}
        for i, v in enumerate(c):
            row[i] = v
        output['detailed'].append(row)
    return output



def write_report(output):
    """
    Write the report
    """
    print("#\n# Total %d calls across %d VCFs\n#" %
          (output['total_calls'], output['num_vcfs']))

    print("#File\tNumCalls")
    for name in output['vcfs']:
        print("%s\t%d" % (name, output['vcf_counts'][name]))

    print("#\n# Summary of consistency\n#")
    print("#VCFs\tCalls\tPct")
    for row in output['shared']:
        print("%d\t%d\t%.2f%%" % (
            row['vcf_count'],
            row['num_calls'],
            row['call_pct'] * 100
        ))

    print("#\n# Breakdown of VCFs' consistency\n#")
    print("#Group\tTotal\tTotalPct\tPctOfFileCalls")
    for row in output['detailed']:
        group = row['group']
        ncalls = row['total']
        tot_pct = row['total_pct'] * 100
        c_text = " ".join(["%.2f%%" % (row[i] * 100)
             if row['group'][i] == "1"
             else "0%"
             for i in range(output['num_vcfs'])
        ])
        print("%s\t%d\t%.2f%%\t%s" %
              (group, ncalls, tot_pct, c_text))

def consistency_main(args):
    """
    Run the program
    """
    args = parse_args(args)

    all_presence, n_calls_per_vcf = read_files(args.allVCFs)
    data = make_report(args.allVCFs, all_presence, n_calls_per_vcf)
    if args.json:
        print(json.dumps(data, indent=4))
    else:
        write_report(data)
