from statsmodels.stats.multitest import fdrcorrection
import numpy as np
import pandas as pd
from scipy import interpolate
from tqdm import tqdm
import output


def testing_procedure(null_dist_scores, input_scores, mode, palpha, ppvalue_correction,
                      poutput, num_known_assoc, known_assoc):

    print("Preparing cdf")
    null_dist_scores_cdf = cdf(null_dist_scores, mode)

    print("Calculating p values")
    p_values_unadj = get_p_values(input_scores, null_dist_scores_cdf)

    print("Determining significant p-values")
    sig_lvl = float(palpha)
    if ppvalue_correction == "none":
        p_values_adj = p_values_unadj
    else:
        alpha_star = sig_lvl
        if num_known_assoc > 0:
            m = count_num_tests(p_values_unadj)
            # m - num_known_assoc is upper bound to m_0
            alpha_star *= m / (m - num_known_assoc)

        if ppvalue_correction == "fdr":
            p_values_adj = fdr(p_values_unadj, alpha_star)

        elif ppvalue_correction == "bonferroni":
            numb_independent_tests = count_num_tests(p_values_unadj)
            sig_lvl = bonferroni_correction(alpha_star, numb_independent_tests)
            p_values_adj = p_values_unadj

    print("Significance level: ", sig_lvl)
    output.write_log(poutput, f"Significance level for {mode}: {sig_lvl}")

    # store p-values of all known associations in case they were not significant
    if known_assoc is not None:
        for (gene1, gene2) in known_assoc:
            try:
                known_assoc[(gene1, gene2)] = (p_values_unadj.loc[gene1, gene2],
                                               p_values_adj.loc[gene1, gene2])
            except KeyError:  # The genes might have been filtered out during preprocessing
                known_assoc[(gene1, gene2)] = (1.0, 1.0)

    # drop all non-associated genes and set p-values > alpha to nan
    p_values_adj, significant_indices = crop_significance(p_values_adj, sig_lvl)

    if significant_indices[0].size == 0:
        exit("No significant gene pairs were found.")

    return significant_indices, p_values_unadj, p_values_adj, sig_lvl, known_assoc


def crop_significance(df, alpha):
    """
    Set p-values > alpha to nan, drop all non-associated genes, and find out indices of significant
    gene pairs

    Parameters
    ----------
    df: pandas dataframe containing gene-gene pairs with their respective p-values
    alpha: significance level

    Returns
    -------
    df : cropped p-value dataframe
    sig_indices : np.array of [row_indices, col_indices], both arrays of indices of significant
                  gene pairs
    """

    # fill every cell > alpha and the diagonal with nan
    df = df[(df < float(alpha))]
    np.fill_diagonal(df.values, np.nan)

    # drop any genes that are not significantly associated with any other gene
    df = df.dropna(axis=0, how="all")
    df = df.dropna(axis=1, how="all")

    # mask will hold True for values in upper triangle
    mask = np.ones(df.shape, dtype='bool')
    mask[np.triu_indices(len(df))] = False

    # get upper triangle indices of significant gene pairs
    sig_indices = df[(df < float(alpha)) & mask].notna().values.nonzero()

    return df, sig_indices


def cdf(scores, assoc_dissoc):
    """
    scores: dictionary containing scores as keys and their respective frequency as values
    assoc_dissoc: str whether testing for "association" or "dissociation" or "both"
    return ecdf array
    """

    # Copy scores dictionary without scores that do not occur
    temp = {}
    for score in scores:
        if scores[score] != 0:
            temp[score] = scores[score]

    # DataFrame with columns 0 (=score) and 1 (=frequency)
    cdf = pd.DataFrame(temp.items())

    # sort by scores
    if assoc_dissoc == "association":
        cdf = cdf.sort_values(by=[0], ascending=False, ignore_index=True)
    else:
        cdf = cdf.sort_values(by=[0], ignore_index=True)

    # Sum up occurrences of all scores and normalize to get frequencies in [0,1]
    cumsum = np.cumsum(list(cdf.iloc[:, 1]))
    totalsum = cumsum[-1]
    cumsum = cumsum/totalsum

    # Interpolate between the scores and insert values beyond the scores that have been simulated
    cdf = interpolate.interp1d(cdf[0], cumsum, bounds_error=False,
                               fill_value=(1, 0) if assoc_dissoc == "association" else (0, 1))

    return cdf


def get_p_values(df, cdf):
    """Calculating p values
    df: pandas dataframe containing gene-gene pairs with their respective scores
    cdf: cumulative distribution function
    return: pandas dataframe containing gene-gene pairs with their respective p values
    """

    temp = np.float32(df.to_numpy())

    # pure numpy operation is faster than pandas map / apply
    for i, x in enumerate(tqdm(temp)):
        temp[i] = cdf(x)

    return pd.DataFrame(data=temp, index=df.index, columns=df.index)


def fdr(df, alpha):
    """false discovery rate
    df: dataframe containing p values of gene gene pairs
    alpha: significance level
    """
    temp = df.copy()
    temp = temp.to_numpy()
    # storing lower diagonal indices
    temp_tril = np.tril_indices(np.size(temp, 0), k=-1)

    # correcting according to yakutiel...
    fdr_values = fdrcorrection(temp[temp_tril], alpha=alpha, method='n')
    temp[temp_tril] = fdr_values[1]

    # highly probable that unnecessary
    # storing upper diagonal indices!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    temp_triu = np.triu_indices(np.size(temp, 0), k=1)
    temp[temp_triu] = temp.T[temp_triu]

    # new dataframe containing only lower triangle p values
    df_fdr = pd.DataFrame(temp, columns=list(df.columns.values), index=list(df.index.values))

    return df_fdr


def bonferroni_correction(alpha, numb_tests):
    """Applying bonferroni correction to significance level
    alpha: significance level
    numb_test: number of independently conducted tests
    return: adjusted significance level
    """
    return alpha/numb_tests


def count_num_tests(df):
    """Determining how many independently conducted tests there are
    df: pandas dataframe of the (currently terminal) scores
    return: number of tests
    """
    n = len(list(df.columns))
    return n * (n-1) / 2
