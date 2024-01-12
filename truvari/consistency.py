"""
Over multiple vcfs, calculate their intersection/consistency.

Calls will match between VCFs if they have a matching key of `CHROM POS ID REF ALT`
"""
import json
import argparse
from collections import defaultdict, Counter

import truvari

def parse_vcf(fn, no_dups=False):
    """
    Simple vcf reader
    """
    fh = truvari.opt_gz_open(fn)
    seen = Counter()
    for line in fh:
        if line.startswith("#"):
            continue
        # Only keep the first 5 fields, and use them as key
        key = "\t".join(line.split("\t")[:5])
        is_dup = key in seen
        if no_dups and is_dup:
            continue
        if is_dup:
            key += f".{str(seen[key])}"
        seen[key] += 1
        yield key


def read_files(allVCFs, no_dups=False):
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
    n_calls_per_vcf = [0] * len(allVCFs)
    for i, vcf in enumerate(allVCFs):
        flag = 1 << (n_vcfs - i - 1)
        for key in parse_vcf(vcf, no_dups):
            # Don't allow duplicates to contribute to the count
            if not all_presence[key] & flag:
                n_calls_per_vcf[i] += 1
                all_presence[key] |= (1 << (n_vcfs - i - 1))

    # We don't care about the calls anyway for stats
    # Then this becomes a list of integers, which will save memory
    return all_presence, n_calls_per_vcf


def parse_args(args):
    """ parse args """
    parser = argparse.ArgumentParser(prog="consistency", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--no-dups", action="store_false",
                        help="Disallow duplicate SVs")
    parser.add_argument("-j", "--json", action='store_true',
                        help="Output report in json format")
    parser.add_argument("-o", "--output", default=None, type=str,
                        help="Write tsv of variant keys and their flag")
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
    output['detailed'] = []

    # setup shared
    for i in reversed(range(n_vcfs)):
        output['shared'].append({'vcf_count': i + 1, 'num_calls': 0, 'call_pct': 0})

    all_overlap = Counter(all_presence)
    for group, ncalls in sorted(all_overlap.items(), key=lambda x: (-x[1], x[0])):
        group = bin(group)
        shared_calls_vcf_count = group.count("1")
        output['shared'][-shared_calls_vcf_count]["num_calls"] += ncalls

        group = group[2:].rjust(n_vcfs, '0')
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

    # finish shared
    for i in range(n_vcfs):
        output['shared'][i]['call_pct'] = output['shared'][i]['num_calls'] / total_unique_calls
    return output



def write_report(output):
    """
    Write the report
    """
    print(f"#\n# Total {output['total_calls']} calls across {output['num_vcfs']} VCFs\n#")

    print("#File\tNumCalls")
    for name in output['vcfs']:
        print(f"{name}\t{output['vcf_counts'][name]}")

    print("#\n# Summary of consistency\n#")
    print("#VCFs\tCalls\tPct")
    for row in output['shared']:
        print(f"{row['vcf_count']}\t{row['num_calls']}\t{row['call_pct'] * 100:.2f}%")

    print("#\n# Breakdown of VCFs' consistency\n#")
    print("#Group\tTotal\tTotalPct\tPctOfFileCalls")
    for row in output['detailed']:
        group = row['group']
        ncalls = row['total']
        tot_pct = row['total_pct'] * 100
        c_text = " ".join([f"{row[i] * 100:.2f}%"
             if row['group'][i] == "1"
             else "0%"
             for i in range(output['num_vcfs'])
        ])
        print(f"{group}\t{ncalls}\t{tot_pct:.2f}%\t{c_text}")

def consistency_main(args):
    """
    Run the program
    """
    args = parse_args(args)

    all_presence, n_calls_per_vcf = read_files(args.allVCFs, args.no_dups)
    data = make_report(args.allVCFs, all_presence.values(), n_calls_per_vcf)
    if args.json:
        for grp in data['detailed']:
            # rename to file
            for idx, name in enumerate(args.allVCFs):
                pct = grp[idx]
                del grp[idx]
                grp[name] = pct
        print(json.dumps(data, indent=4))
    else:
        write_report(data)
    if args.output:
        with open(args.output, 'w') as fout:
            for k,v in all_presence.items():
                fout.write(f"{k}\t{v}\n")
