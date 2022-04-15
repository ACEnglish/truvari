"""
Annotate low complexity region entropy score for variants
"""
import re
import sys
import math
import pysam

"""
Credit
https://jszym.com/blog/dna_protein_complexity/
"""

def overlapping_windows(sequence, L):
    """
    Returns overlapping windows of size `L` from sequence `sequence`
    :param sequence: the nucleotide or protein sequence to scan over
    :param L: the length of the windows to yield
    """
    windows = []

    for index, residue in enumerate(sequence):
        if (index + L) < (len(sequence) + 1):
            window = sequence[index:L+index]
            windows.append(window)

    return windows

def sequence_to_repvec(sequence, N):
    """
    Computes the repetition vector (as seen in Wooton, 1993) from a
    given sequence of a biopolymer with `N` possible residues.

    :param sequence: the nucleotide or protein sequence to generate a repetition vector for.
    :param N: the total number of possible residues in the biopolymer `sequence` belongs to.
    """
    encountered_residues = set()
    repvec = []

    for residue in sequence:
        if residue not in encountered_residues:
            residue_count = sequence.count(residue)

            repvec.append(residue_count)

            encountered_residues.add(residue)

        if len(encountered_residues) == N:
            break

    while len(repvec) < N:
        repvec.append(0)

    return sorted(repvec, reverse=True)

def sequence_entropy(sequence, N=4):
    """
    Computes the Shannon Entropy of a given sequence of a
    biopolymer with `N` possible residues. See (Wooton, 1993)
    for more.

    :param sequence: the nucleotide or protein sequence whose Shannon Entropy is to calculated.
    :param N: the total number of possible residues in the biopolymer `sequence` belongs to.
    """
    repvec = sequence_to_repvec(sequence, N)

    L = len(sequence)

    entropy = sum([-1*(n/L)*math.log((n/L), N) for n in repvec if n != 0])

    return entropy

def mask_low_complexity(seq, L=12, N=20, maskchar="x"):
    """
    Mask LCR sequence
    """
    windows = overlapping_windows(seq, L)

    rep_vectors = [(window, compute_rep_vector(window, N)) for window in windows]

    window_complexity_pairs = [(rep_vector[0], complexity(rep_vector[1], N)) for rep_vector in rep_vectors]

    complexities = np.array([complexity(rep_vector[1], N) for rep_vector in rep_vectors])

    avg_complexity = complexities.mean()
    std_complexity = complexities.std()

    k1_cutoff = min([avg_complexity + std_complexity,
                 avg_complexity - std_complexity])

    alignment = [[] for i in range(0, len(seq))]

    for window_offset, window_complexity_pair in enumerate(window_complexity_pairs):

        if window_complexity_pair[1] < k1_cutoff:
            window = "".join([maskchar for i in range(0, L)])

        else:
            window = window_complexity_pair[0]

        for residue_offset, residue in enumerate(window):
            i = window_offset+residue_offset
            alignment[i].append(residue)

    new_seq = []

    for residue_array in alignment:
        if residue_array.count(maskchar) > 3:
            new_seq.append(maskchar)
        else:
            new_seq.append(residue_array[0])

    new_seq = "".join(new_seq)

    return (new_seq, alignment)

 def edit_header(my_vcf):
    """
    Add INFO for new field to vcf
    """
    header = my_vcf.header.copy()
    desc = 'Description="' + ",".join([f"{_}x" for _ in bins]) + '">'
    header.add_line(('##INFO=<ID=LCR,Number=1,Type=Float,'
                     'Description="Low complexity region entropy score">'))
    return header

def add_lcr(vcf, n_header):
    """
    Adds LCR to each entry in VCF and yields them
    """
    if not n_header:
        n_header = edit_header(vcf)
    for entry in vcf:
        seq = entry.alts[0] if len(entry.alts[0]) > len(entry.ref) else entry.ref
        score = sequence_entropy(seq)
        entry.translate(n_header)
        entry.info["LCR"] = score
        yield entry

def lcr_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)
    vcf = pysam.VariantFile(args.input)
    n_header = edit_header(vcf, bins)
    out = pysam.VariantFile(args.output, 'w', header=n_header)
    for entry in add_lcr(vcf, n_header):
        out.write(entry)
    out.close()
    logging.info("Finished lcr")
