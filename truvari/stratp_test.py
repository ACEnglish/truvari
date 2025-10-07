"""
Stratification Performance Test
"""
import os
import sys
import logging
import argparse
import itertools

import joblib
import numpy as np
import pandas as pd
import truvari
from truvari.vcf2df import vcf2df_main


def fdr_bh(pvals, alpha=0.05):
    """
    Benjamini-Hochberg FDR correction. Replicates statsmodels.stats.multitest bh

    lifted from
    Parameters:
    pvals : array-like
        List or array of p-values.
    alpha : float
        Desired FDR control level.

    Returns:
    rejected : np.ndarray
        Boolean array where True means the hypothesis is rejected.
    pvals_corrected : np.ndarray
        Adjusted p-values.
    """
    pvals = np.asarray(pvals)
    n = len(pvals)
    sorted_indices = np.argsort(pvals)
    sorted_pvals = pvals[sorted_indices]

    adjusted = sorted_pvals * n / (np.arange(1, n + 1))
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    rejected = adjusted <= alpha

    # Return results in original order
    pvals_corrected = np.empty_like(adjusted)
    pvals_corrected[sorted_indices] = adjusted
    rejected_final = np.empty_like(rejected)
    rejected_final[sorted_indices] = rejected

    return rejected_final, pvals_corrected


def check_co_occur(data, features, threshold=0.50):
    """
    Calculate phi to detemine if there's a strong co-occurance of feature values
    If abs(phi) >= 0.50, the features could be redundant for stratp test
    """
    # pylint: disable=import-outside-toplevel
    try:
        from scipy.stats import chi2_contingency
    except (OSError,  ModuleNotFoundError):
        logging.error("Cannot run `--check-co` without scipy package")
        sys.exit(1)
    for a, b in itertools.combinations(features, 2):
        if data[a].nunique() < data[b].nunique():
            a, b = b, a
        grp = data.groupby([a, b]).size().unstack()
        chi2, _, _, _ = chi2_contingency(grp)
        n = grp.sum().sum()
        phi = np.sqrt(chi2 / n)
        logging.debug("Co-occurance check %s - %s phi=%.3f\n%s",
                      a, b,  phi, grp)
        if abs(phi) >= threshold:
            logging.warning(
                "Features %s and %s may be redundant (phi=%.3f)", a, b, phi)


def get_scores(df, state, features, min_obs=10):
    """
    Calculate enrichment scores
    """
    tot_counts = []
    for feat in features:
        cnt = df[feat].value_counts()
        cnt.name = "count"
        cnt = pd.DataFrame(cnt)
        cnt['feat'] = feat
        tot_counts.append(cnt)
    tot_counts = pd.concat(tot_counts)
    tot_counts = tot_counts.reset_index().set_index(['feat', 'index'])

    logging.debug("Full Counts\n%s", str(tot_counts))
    all_cont_obs = pd.crosstab(df[state], [df[i] for i in features])
    logging.info("%d stratifications had ≥1 observed variant",
                 all_cont_obs.shape[1])
    filtered = all_cont_obs.sum(axis=0) < 10
    if filtered.any():
        logging.warning(
            "%d stratifications with < %d observations ignored", filtered.sum(), min_obs)
        with pd.option_context('display.max_rows', 100):
            logging.debug("Fitered Bins\n%s", str(all_cont_obs.T[filtered]))

    cont_obs = all_cont_obs.T[~filtered].T
    # Calculate Expected
    c_total = cont_obs.sum(axis=0)
    r_total = cont_obs.sum(axis=1)
    g_total = c_total.sum()
    exp = np.outer(r_total, c_total) / g_total

    score = (np.sign(cont_obs.values[1] - exp[1])) * \
            (cont_obs.values[1] - exp[1])**2 / exp[1]

    df_rank = pd.DataFrame(list(zip(list(cont_obs.columns),
                                    list(cont_obs.values[1]), list(exp[1]), list(score))),
                           columns=["Values", "observed", "expected", "score"])

    df_rank.sort_values(["score"], ascending=True, inplace=True)
    df_rank.reset_index(inplace=True, drop=True)
    df_rank.index.names = ['rank']
    df_rank[features] = pd.DataFrame(df_rank["Values"].tolist(),
                                     index=df_rank.index)
    df_rank = df_rank.drop("Values", axis=1)
    return df_rank, tot_counts, all_cont_obs


def permutation_test(a_values, b_values, n_samps=10000, tailed="two"):
    """
    a_values - b_values mean difference from permutations' mean differences
    tail can be "two", "left", or "right"
    """
    delta = a_values.mean() - b_values.mean()
    if tailed == "two":
        delta = abs(delta)
    logging.debug("Case %.2f (%d)", a_values.mean(), len(a_values))
    logging.debug("Control %.2f (%d)", b_values.mean(), len(b_values))
    logging.debug("delta %.2f", delta)

    def run_t(pooled, sizeA, sizeB):
        np.random.shuffle(pooled)
        starA = pooled[:sizeA]
        starB = pooled[-sizeB:]
        return starA.mean() - starB.mean()
    pooled = np.hstack([a_values, b_values])
    estimates = np.array(
        [run_t(pooled, a_values.size, b_values.size) for _ in range(n_samps)])
    if tailed == "left":
        diffCount = len(np.where(estimates <= delta)[0])
    elif tailed == "right":
        diffCount = len(np.where(estimates >= delta)[0])
    else:
        estimates = np.abs(estimates)
        diffCount = len(np.where(estimates >= delta)[0])
    pct_extreme = float(diffCount) / float(n_samps)
    one, fif, nine = np.quantile(estimates, [0.01, 0.5, 0.99])
    logging.debug("Quantiles: 1%% %.4f - 50%% %.4f - 99%% %.4f", one, fif, nine)
    logging.debug("Pct at least as extreme %.2f", pct_extreme)
    return delta, one, fif, nine, pct_extreme


def stratp_test(m_data, state, features,  min_obs=10, n_perm=10000, tail="left", alpha=0.01):
    """
    Create the StratP score data-frame
    """
    # pylint: disable=unsubscriptable-object
    rank, counts, all_counts = get_scores(m_data, state, features, min_obs)
    rows = []

    for feat in features:
        for val in rank[feat].unique():
            sub = rank[feat] == val
            if sub.sum() < min_obs or (~sub).sum() < min_obs:
                logging.warning("Unable to permute %s : %s with < %d groups (%d case, %d control)",
                                feat, val, min_obs, sub.sum(), (~sub).sum())
                result = [np.array(rank[sub].index).mean() - np.array(rank[~sub].index).mean(),
                          np.nan, np.nan, np.nan, np.nan]
            else:
                logging.debug("Testing %s : %s", feat, val)
                result = permutation_test(np.array(rank[sub].index),
                                        np.array(rank[~sub].index),
                                        n_perm, tail)
            acc = m_data[m_data[feat] == val][state].mean()
            rows.append([feat,
                         val,
                         counts.loc[feat, val]["count"],
                         acc,
                         *result])

    plt_rank = pd.DataFrame(rows, columns=["feature", "value", "obs",
                                           "acc", "delta", "q1", "q50",
                                           "q99", "pval"])
    reject, adj_pval = fdr_bh(pvals=plt_rank['pval'], alpha=alpha)
    adj = pd.concat([plt_rank,
                     pd.Series(adj_pval, name="adj_pval"),
                     pd.Series(reject, name="reject"),
                     ], axis=1)
    # pylint: enable=unsubscriptable-object
    return adj, all_counts


def parse_args(args):
    """
    Parse args
    """
    parser = argparse.ArgumentParser(prog="stratp", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bench_dir", metavar="DIR", type=str,
                        help="Truvari bench result directory")
    parser.add_argument("--features", type=str, default=None,
                        help="Comma-separated list of categorical features to score (e.g. svtype,szbin)")
    # TODO: could do f1 as the accuracy score for both
    parser.add_argument("--states", type=str, default='base', choices=['base', 'comp'],
                        help="Score VCF entries from (base)line or (comp)arison (default:base)")
    parser.add_argument("--preset", type=str, default=None, choices=PRESETS.keys(),
                        help="Preset features to score")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Score output tsv file (stdout)")
    parser.add_argument("-c", "--count-output", type=str, default=None,
                        help="Stratification variant count tsv file")
    parser.add_argument("--debug", action='store_true',
                        help="Verbose logging")

    thresg = parser.add_argument_group("statistics arguments")
    thresg.add_argument("--min-obs", type=int, default=10,
                        help="Minimum observation for a bin (%(default)s)")
    thresg.add_argument("--nperm", type=int, default=10000,
                        help="Number of permutations to perform (%(default)s)")
    thresg.add_argument("--tail", type=str, default="left", choices=['left', 'right', 'two'],
                        help="Tail to test (%(default)s)")
    thresg.add_argument("--alpha", type=truvari.restricted_float, default=0.01,
                        help="Significance threshold (%(default)s)")
    thresg.add_argument("--check-co", type=float, default=-1,
                        help="Check if features co-occur and warns if phi ≥ threshold. Requires scipy. (-1=off)")
    thresg.add_argument("--seed", type=int, default=42,
                        help="Random seed for permutation tests (%(default)s)")

    args = parser.parse_args(args)
    truvari.setup_logging(args.debug, show_version=False)
    if (args.features is None) == (args.preset is None):
        logging.error("Either `--features` xor `--preset` must be provided")
        sys.exit(1)
    return args


def set_mims_bins(data):
    """
    Edits the DataFrame in place to setup MIMS SV benchmark stratification
    Then returns the feature column names and state column to use in the scoring
    """
    data['is_tp'] = data['state'].str.startswith('tp')
    data['isolated'] = data['NumNeighbors'] == 0
    bins = [0, 0.05, 0.3, 0.7, 1]
    labels = ["SomaticLow(<5%)", "Somatic", "GermlineHet", "GermlineHom"]
    data['VAF_bin'] = pd.cut(data['VAF_alt'], bins=bins, labels=labels)
    return ['svtype', 'szbin', 'TRF', 'isolated', 'VAF_bin'], 'base'


def set_giabv1_bins(data):
    """
    Edits the DataFrame in place to setup GIAB v1.1 SV benchmark stratification
    Then returns the feature column names and state column to use in the scoring
    """
    data['is_tp'] = data['state'].str.startswith('tp')
    data['is_inter'] = ~data['RM_clsfam'].isna()
    data['GT'] = data['HG002_GT'].apply(lambda x: truvari.get_gt(x).name)
    return ['svtype', 'szbin', 'TRF', 'is_inter', 'GT'], 'base'


def set_basic_bins(data):
    """
    Edits the DataFrame in place to setup basic SV benchmark stratification
    Then returns the feature column names and state column to use in the scoring
    TODO
    This doesn't exactly work because there's too few groups to do permutation testing on the szbin
    e.g. szbin : [50,100) with < 5 observations ignored (4 case, 36 control)
    I could wrap up lcr and annotate the alleles on the fly because that's essentially tandem repeats
    But then I'd need to work on the vcf2df_main to get the alleles in the DataFrame
    And I'm not positive that wouldn't be slow
    But, also, if you do lcr, why not do NumNeighbors, also, and you think, well, we could ask users to run those before
    hand, but then why not ask them to run all the annotations. And then why not automate all of that work in the tool
    """
    data['is_tp'] = data['state'].str.startswith('tp')
    first_gt = [_ for _ in data.columns if _.endswith('GT')][0]
    data['GT'] = data[first_gt].apply(lambda x: truvari.get_gt(x).name)
    return ['svtype', 'szbin', 'GT'], 'base'


# Presets
PRESETS = {
    'mims': set_mims_bins,
    'giab1.1': set_giabv1_bins,
    # 'basic': set_basic_bins,
}


def stratp_main(args):
    """
    Run StratP test on a truvari directory with benchmarking results
    """
    args = parse_args(args)
    np.random.seed(args.seed)
    data_file = os.path.join(args.bench_dir, "data.jl")
    if not os.path.exists(data_file):
        logging.info("Constructing DataFrame in %s", data_file)
        vcf2df_main(f"-i -f -b {args.bench_dir} {data_file}".split(' '))

    data = joblib.load(data_file)
    if args.preset:
        args.features, args.states = PRESETS[args.preset](data)
    else:
        args.features = args.features.split(',')
        data['is_tp'] = data['state'].str.startswith('tp')

    if args.states == 'base':
        data = data[data['state'].isin(['tpbase', 'fn'])]
    elif args.states == "comp":
        data = data[data['state'].isin(['tp', 'fp'])]

    # Bin syntax NAME:[f|w]
    for idx, col in enumerate(args.features):
        if ':' in col:
            name, key = col.split(':')
            method, n_bins = key[0], key[1:]
            if method == 'f':
                n_bins = int(n_bins)
                data[name] = pd.qcut(data[name], q=n_bins, duplicates='drop')
            elif method == 'w':
                n_bins = int(n_bins)
                data[name] = pd.cut(data[name], bins=n_bins,
                                    right=False, include_lowest=True)
            elif method == 'c':
                try:
                    bins = list(map(int, n_bins.split('-')))
                except ValueError:
                    bins = list(map(float, n_bins.split('-')))
                data[name] = pd.cut(data[name], bins=bins,
                                    right=False, include_lowest=True)
            args.features[idx] = name

    n_strats = 1
    n_values = 0
    for i in args.features:
        data[i] = data[i].astype(str)
        c = data[i].nunique()
        n_strats *= c
        n_values += c
    n_feat = len(args.features)

    logging.info("Scoring %d features with %d values (max %d stratifications)",
                 n_feat, n_values, n_strats)

    if args.check_co != -1:
        check_co_occur(data, args.features, args.check_co)

    scores, counts = stratp_test(data, state="is_tp",
                                 features=args.features,
                                 min_obs=args.min_obs,
                                 n_perm=args.nperm,
                                 tail=args.tail,
                                 alpha=args.alpha)

    if args.count_output:
        counts.T.to_csv(args.count_output, sep='\t')
    scores.round(4).sort_values(by=['adj_pval', 'pval', 'acc'], ascending=True).to_csv(
        args.output, sep='\t', index=False)
    logging.info("Finished")


if __name__ == '__main__':
    stratp_main(sys.argv[1:])
