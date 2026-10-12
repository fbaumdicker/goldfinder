import numpy as np
from tqdm import tqdm
import pandas as pd
from numba import njit, prange
from numba.typed import Dict
from numba import types
from numba_progress import ProgressBar
import output


def scoring_procedure(pscore, simul_anc_states, simul_desc_states, simul_leaves, simul_counts,
                      fitch_score_root_one, fitch_score_root_zero, fitch_score_root_both,
                      tree_struc, branch_distances, pgenes_simulated, pre_df, input_anc_state,
                      input_desc_state, input_leaves_state, poutput):
    """
    pscore: string describing which score to use
    """

    if pscore == "simultaneous":

        print("Calculating simultaneous score for the null distribution")
        with ProgressBar(total=len(simul_anc_states)) as progress:
            null_dist_scores_numba = simultaneous_score_null_dist_numba(
                np.asarray(simul_anc_states), np.asarray(simul_desc_states),
                np.asarray(simul_counts), progress)

        # transform scores array into normal python dictionary
        null_dist_scores = {}
        num_ancs = np.asarray(simul_anc_states).shape[1]
        for i, counts in enumerate(null_dist_scores_numba):
            # see comments in simultaneous_score_null_dist_numba on what the index encodes
            null_dist_scores[i - num_ancs] = counts

        output.write_distribution(null_dist_scores, poutput + "/" +
                                  "score_distribution_before_additional_simulation.txt")

        print("Calculating simultaneous score for the input data")
        with ProgressBar(total=len(input_anc_state)) as progress:
            input_scores = simultaneous_score_input_numba(np.asarray(input_anc_state),
                                                          np.asarray(input_desc_state), progress)
        input_scores = pd.DataFrame(input_scores, columns=list(pre_df.index),
                                    index=list(pre_df.index))

    if pscore == "terminal":
        null_dist_scores = terminal_score_null_dist(simul_counts, simul_leaves)
        input_scores = terminal_score_input(pre_df, input_leaves_state)

    if pscore == "subsequent":
        null_dist_scores = subsequent_score_null_dist(simul_counts, simul_anc_states,
                                                      simul_desc_states)
        input_scores = subsequent_score_input(pre_df, input_anc_state, input_desc_state)

    if pscore == "coinfinder":
        null_dist_scores = coinfinder_score_null_dist(simul_counts, simul_leaves)
        input_scores = coinfinder_score_input(pre_df, input_leaves_state)

    return null_dist_scores, input_scores


def additional_scoring(add_simul_anc_states, add_simul_desc_states, add_simul_counts, k, k_frac,
                       null_dist_scores, poutput):
    """


    Parameters
    ----------
    add_simul_anc_states : TYPE
        DESCRIPTION.
    add_simul_desc_states : TYPE
        DESCRIPTION.
    add_simul_counts : TYPE
        DESCRIPTION.
    k : int
        score threshold of 'extreme' scores
    k_frac : float
        fraction of Fitch scores >= k

    Returns
    -------
    None.

    """
    # Calculate scores of additional simulated states

    print("Calculating simultaneous score for the new null distribution")
    with ProgressBar(total=len(add_simul_anc_states)) as progress:
        add_null_dist_scores_numba = simultaneous_score_null_dist_numba(
            np.asarray(add_simul_anc_states), np.asarray(add_simul_desc_states),
            np.asarray(add_simul_counts),  progress)

    # transform scores array into normal python dictionary
    add_null_dist_scores = {}
    num_ancs = np.asarray(add_simul_anc_states).shape[1]
    for i, counts in enumerate(add_null_dist_scores_numba):
        add_null_dist_scores[i - num_ancs] = counts

    # Scaling the results to pretend we simulated over the whole score range.
    add_null_dist_scores = {key: add_null_dist_scores[key] * k_frac ** 2 for key in add_null_dist_scores}

    # Add values of null_dist_scores for keys < k and more precise values of add_null_dist_scores
    # for keys >= k and keys <= -k into an updated dictionary
    updated_null_dist_scores = {}
    for key in null_dist_scores:
        if null_dist_scores[key] != 0 and key < k:
            updated_null_dist_scores[key] = null_dist_scores[key]

    for key in add_null_dist_scores:
        if add_null_dist_scores[key] != 0 and (key >= k or key <= -k):
            updated_null_dist_scores[key] = add_null_dist_scores[key]

    output.write_distribution(updated_null_dist_scores,
                              f"{poutput}/score_distribution_after_additional_simulation.txt")

    return updated_null_dist_scores

###################################################################################################


def coinfinder_score(g, p):
    score = np.sum(np.multiply(g, p), axis=1)/len(g)
    return score


def coinfinder_score_null_dist(counts, leaf_states):
    """Determining the coinfinder scores for the null distribution
    counts: list containing how often a specific pattern of in leaf_states occurred in the simulated
            data
    leaf_states: list containing the gene presence/absence state at the leaf nodes
    return: dictionary with score as key and its respective frequency as value
    """
    print("Calculating coinfinder score for the null distribution")
    counts = np.asarray(counts)
    leaves = np.asarray(leaf_states)
    coinfinder_scores = {}
    for x in tqdm(range(len(leaves))):
        c_score = coinfinder_score(leaves[x], leaves)
        new_counts = counts * counts[x]
        new_counts[x] = counts[x] * (counts[x] - 1) / 2
        for y in range(len(c_score)):
            if y >= x:
                coinfinder_scores[c_score[y]] = coinfinder_scores.get(c_score[y], 0) + new_counts[y]
    return coinfinder_scores


def coinfinder_score_input(df, leaf_states):
    """Determining the coinfinder scores for the input data
    df: input gene presence/absence data
    leaf_states: list containing the gene presence/absence state at the leaf nodes
    return: pandas data frame of gene-gene pairs with their respective terminal score
    """
    print("Calculating coinfinder scores for the input data")
    leaves = np.asarray(leaf_states)
    gene_names = list(df.index)
    df_scores = pd.DataFrame(columns=gene_names, index=gene_names)

    for x in tqdm(range(len(leaves))):
        c_score = coinfinder_score(leaves[x], leaves)
        df_scores.iloc[x] = c_score

    return df_scores

####################################################################################################


def terminal_score(g, p):
    """Measuring sample-wide association across leaves of phylogenetic tree
    g: gene presence/absence state at leaves of one gene as np.array
    p: gene presence/absence state at leaves of all genes as np.array
    return: np.array of scores of one gene compared against all other genes
    """
    score = np.sum(np.multiply(g, p)*4 + 1 + np.negative(p)*2 + np.negative(g)*2, axis=1) / len(g)
    return score


def terminal_score_null_dist(counts, leaf_states):
    """Determining the terminal scores for the null distribution
    counts: list containing how often a specific pattern of in leaf_states occurred in the simulated
            data
    leaf_states: list containing the gene presence/absence state at the leaf nodes
    return: dictionary with score as key and its respective frequency as value
    """
    print("Calculating terminal score for the null distribution")
    counts = np.asarray(counts)
    leaves = np.asarray(leaf_states)
    terminal_scores = {}
    for x in tqdm(range(len(leaves))):
        t_score = terminal_score(leaves[x], leaves)
        new_counts = counts * counts[x]
        new_counts[x] = counts[x] * (counts[x] - 1) / 2
        for y in range(len(t_score)):
            if y >= x:
                terminal_scores[t_score[y]] = terminal_scores.get(t_score[y], 0) + new_counts[y]
    return terminal_scores


def terminal_score_two_null_dists(counts_one, leaf_states_one, counts_two, leaf_states_two):
    """Determining the terminal scores for the null distribution
    counts: list containing how often a specific pattern of in leaf_states occurred in the simulated
            data
    leaf_states: list containing the gene presence/absence state at the leaf nodes
    return: dictionary with score as key and its respective frequency as value
    """
    print("Calculating terminal score for the null distribution")
    counts_one = np.asarray(counts_one)
    leaves_one = np.asarray(leaf_states_one)

    counts_two = np.asarray(counts_two)
    leaves_two = np.asarray(leaf_states_two)

    terminal_scores = {}
    for x in tqdm(range(len(leaves_one))):
        t_score = terminal_score(leaves_one[x], leaves_two)
        new_counts = counts_two * counts_one[x]

        for y in range(len(t_score)):
            # if y >= x:
            terminal_scores[t_score[y]] = terminal_scores.get(t_score[y], 0) + new_counts[y]
    return terminal_scores


def terminal_score_input(df, leaf_states):
    """Determining the terminal scores for the input data
    df: input gene presence/absence data
    leaf_states: list containing the gene presence/absence state at the leaf nodes
    return: pandas data frame of gene-gene pairs with their respective terminal score
    """
    print("Calculating terminal scores for the input data")
    leaves = np.asarray(leaf_states)
    gene_names = list(df.index)
    df_scores = pd.DataFrame(columns=gene_names, index=gene_names)

    for x in tqdm(range(len(leaves))):
        t_score = terminal_score(leaves[x], leaves)
        df_scores.iloc[x] = t_score

    return df_scores


def temporay_terminal_score(df):
    terminal_scores = {}
    num_genomes = len(df.columns)
    df = df.to_numpy()
    for x in tqdm(range(len(df))):
        terminal_score = np.sum(np.multiply(
            df[x], df)*4 + 1 + np.negative(df)*2 + np.negative(df[x])*2, axis=1) / num_genomes
        for y in range(len(terminal_score)):
            if y >= x:
                terminal_scores[terminal_score[y]] = terminal_scores.get(terminal_score[y], 0) + 1
    return terminal_scores

####################################################################################################


def simultaneous_score(panc, pdesc, ganc, gdes):
    """Measuring degree of parallel change in one gene against another gene across branches of
       phylogenetic tree
    g: gene presence/absence state at leaves of one gene as np.array
    p: gene presence/absence state at leaves of all genes as np.array
    return: np.array of scores of one gene compared against all other genes
    """
    score = np.sum(np.multiply((panc-pdesc), (ganc-gdes)), axis=1)
    return score


def simultaneous_score_null_dist(anc, desc, counts):
    """Determining the simultaneous scores for the null distribution
    counts: list containing how often a specific pattern of in leaf_states occurred in the simulated
            data
    anc: list containing the gene presence/absence state at the ancestral node of the desc node
    desc: list containing the gene presence/absence state at a node
    return: dictionary with score as key and its respective frequency as value
    """

    print("Calculating simultaneous score for the null distribution")
    counts = np.asarray(counts)
    ancestors = np.asarray(anc)
    descendants = np.asarray(desc)
    simu_scores = {}

    for x in tqdm(range(len(ancestors))):
        simu_score = simultaneous_score(
            ancestors[x:], descendants[x:], ancestors[x], descendants[x])
        new_counts = counts[x:] * counts[x]
        new_counts[0] = counts[x] * (counts[x] - 1) / 2
        for y in range(len(simu_score)):
            simu_scores[simu_score[y]] = simu_scores.get(simu_score[y], 0) + new_counts[y]

    return simu_scores


@njit(parallel=True, fastmath=True)
def simultaneous_score_null_dist_numba(anc, desc, counts, progress_proxy):
    # The simultaneous score (as defined by treeWAS) is implemented as a comparison between two genotypes,
    # instead of between a genotype and a phenotype

    # 1D array of number of occurrences of the different patterns
    counts = np.asarray(counts)
    ancestors = np.asarray(anc)
    descendants = np.asarray(desc)

    # ancestors.shape[1] is number of nodes-1 because root has no ancestor
    num_ancs = ancestors.shape[1]

    # simultaneous score can be num_ancs at max, negative this number at min, and also 0
    # use matrix so different threads won't interfere with each other
    score_occurrences = np.zeros((ancestors.shape[0], num_ancs * 2 + 1))

    def idx(curr_pattern_idx, score):
        return (curr_pattern_idx, score + num_ancs)

    for x in prange(len(ancestors)):

        # Performs score calculations of genotype x with all genotypes x, x+1, ..., n at once
        curr_scores = np.sum(np.multiply((ancestors[x:]-descendants[x:]),
                                         (ancestors[x]-descendants[x])), axis=1)

        # each genotype is actually a pattern that might have occured multiple times
        # new_counts is an array of how often the pairs would have been observed if we didnt
        # summarize in patterns
        new_counts = counts[x:] * counts[x]

        # genotypes of pattern x would have been observed n*(n-1)/2 times pairing with another
        # genotype of pattern x
        new_counts[0] = counts[x] * (counts[x] - 1) / 2
        for score_i in range(len(curr_scores)):
            score_occurrences[idx(x, curr_scores[score_i])] += new_counts[score_i]

        progress_proxy.update(1)
    return score_occurrences.sum(axis=0)  # sum results of parallel for loop


@njit(parallel=True, fastmath=True)
def simultaneous_score_null_dist_numba_dict(anc, desc, counts, progress_proxy):
    # Old implementation with dicts. This sometimes threw KeyError, but might be helpful in the
    # future
    # The simultaneous score (as defined by treeWAS) is implemented as a comparison between two genotypes,
    # instead of between a genotype and a phenotype

    # 1D array of number of occurrences of the different patterns
    counts = np.asarray(counts)
    ancestors = np.asarray(anc)
    descendants = np.asarray(desc)

    # dict of score: number of occurrences
    simu_score_dict = Dict.empty(key_type=types.int64, value_type=types.int64)

    for x in prange(len(ancestors)):
        # Performs score calculations of genotype x with all genotypes x, x+1, ..., n at once
        curr_scores = np.sum(np.multiply((ancestors[x:]-descendants[x:]),
                                         (ancestors[x]-descendants[x])), axis=1)

        # each genotype is actually a pattern that might have occured multiple times
        # new_counts is an array of how often the pairs would have been observed if we didnt
        # summarize in patterns
        new_counts = counts[x:] * counts[x]

        # genotypes of pattern x would have been observed n*(n-1)/2 times pairing with another
        # genotype of pattern x
        new_counts[0] = counts[x] * (counts[x] - 1) / 2
        for y in range(len(curr_scores)):
            simu_score_dict[curr_scores[y]] = simu_score_dict.get(curr_scores[y], 0) + new_counts[y]
        progress_proxy.update(1)
    return simu_score_dict


@njit(parallel=True, fastmath=True)
def simultaneous_score_input_numba(anc, desc, progress_proxy):
    ancestors = np.asarray(anc)
    descendants = np.asarray(desc)
    df_scores = np.zeros(shape=(len(anc), len(anc)))

    for x in prange(len(ancestors)):
        df_scores[x] = np.sum(np.multiply((ancestors-descendants),
                              (ancestors[x]-descendants[x])), axis=1)
        progress_proxy.update(1)

    return df_scores


def simultaneous_score_input(df, anc, desc):
    """Determining the simultaneous scores for the input data
    df: input gene presence/absence data
    anc: list containing the gene presence/absence state at the ancestral node of the desc node
    desc: list containing the gene presence/absence state at a node
    return: pandas data frame of gene-gene pairs with their respective simultaneous score
    """
    print("Calculating simultaneous scores for the input data")
    ancestors = np.asarray(anc)
    descendants = np.asarray(desc)
    gene_names = list(df.index)
    df_scores = pd.DataFrame(columns=gene_names, index=gene_names)

    for x in tqdm(range(len(ancestors))):
        simu_score = simultaneous_score(ancestors, descendants, ancestors[x], descendants[x])
        df_scores.iloc[x] = simu_score

    return df_scores


"""
@njit(parallel=True, fastmath=True)
def simultaneous_score_input_numba(anc, desc, progress_proxy):
    Determining the simultaneous scores for the input data
    df: input gene presence/absence data
    anc: list containing the gene presence/absence state at the ancestral node of the desc node
    desc: list containing the gene presence/absence state at a node
    return: pandas data frame of gene-gene pairs with their respective simultaneous score

    print("Calculating simultaneous scores for the input data")
    ancestors = np.asarray(anc)
    descendants = np.asarray(desc)
    #gene_names = list(df.index)
    #df_scores = pd.DataFrame(columns = gene_names, index = gene_names)
    df_scores = np.empty((len(ancestors), len(ancestors)))


    for x in prange(len(ancestors)):
        simu_score = np.sum(np.multiply((ancestors-descendants),(ancestors[x]-descendants[x])),
                            axis=1)
        #simultaneous_score(ancestors, descendants, ancestors[x],descendants[x])
        #df_scores[x] = simu_score
        df_scores[x] = simu_score
        progress_proxy.update(1)
    print(df_scores)
    return df_scores
"""

####################################################################################################


def subsequent_score(panc, pdesc, ganc, gdes):
    """Measuring proportion of tree in which two genes co-exist
    g: gene presence/absence state at leaves of one gene as np.array
    p: gene presence/absence state at leaves of all genes as np.array
    return: np.array of scores of one gene compared against all other genes
    """
    score = np.sum(np.multiply(panc, ganc)*4/3 + np.multiply(panc, gdes)*2/3 + np.multiply(pdesc,
                   ganc)*2/3 + np.multiply(pdesc, gdes)*4/3 - panc-pdesc-ganc-gdes + 1, axis=1)
    return score


def subsequent_score_null_dist(counts, anc, desc):
    """Determining the subsequent scores for the null distribution
    counts: list containing how often a specific pattern of in leaf_states occurred in the simulated
            data
    anc: list containing the gene presence/absence state at the ancestral node of the desc node
    desc: list containing the gene presence/absence state at a node
    return: dictionary with score as key and its respective frequency as value
    """
    print("Calculating subsequent score for the null distribution")
    counts = np.asarray(counts)
    ancestors = np.asarray(anc)
    descendants = np.asarray(desc)
    sub_scores = {}

    for x in tqdm(range(len(ancestors))):
        sub_score = subsequent_score(ancestors, descendants, ancestors[x], descendants[x])
        new_counts = counts * counts[x]
        new_counts[x] = counts[x] * (counts[x] - 1) / 2
        for y in range(len(sub_score)):
            if y >= x:
                sub_scores[sub_score[y]] = sub_scores.get(sub_score[y], 0) + new_counts[y]
    return sub_scores


def sub_score_two_null_dists(counts_one, anc_one, desc_one, counts_two, anc_two, desc_two):
    """Determining the terminal scores for the null distribution
    counts: list containing how often a specific pattern of in leaf_states occurred in the simulated
            data
    leaf_states: list containing the gene presence/absence state at the leaf nodes
    return: dictionary with score as key and its respective frequency as value
    """
    print("Calculating subsequent score for the null distribution")
    counts_one = np.asarray(counts_one)
    anc_one = np.asarray(anc_one)
    desc_one = np.asarray(desc_one)

    counts_two = np.asarray(counts_two)
    anc_two = np.asarray(anc_two)
    desc_two = np.asarray(desc_two)

    sub_scores = {}
    for x in tqdm(range(len(anc_one))):
        s_score = subsequent_score(anc_one[x], desc_one[x], anc_two, desc_two)
        new_counts = counts_two * counts_one[x]

        for y in range(len(s_score)):
            # if y >= x:
            sub_scores[s_score[y]] = sub_scores.get(s_score[y], 0) + new_counts[y]
    return sub_scores


def subsequent_score_input(df, anc, desc):
    """Determining the simultaneous scores for the input data
    df: input gene presence/absence data
    anc: list containing the gene presence/absence state at the ancestral node of the desc node
    desc: list containing the gene presence/absence state at a node
    return: pandas data frame of gene-gene pairs with their respective subsequent score
    """
    print("Calculating subsequent scores for the input data")
    ancestors = np.asarray(anc)
    descendants = np.asarray(desc)
    gene_names = list(df.index)
    df_scores = pd.DataFrame(columns=gene_names, index=gene_names)

    for x in tqdm(range(len(ancestors))):
        sub_score = subsequent_score(ancestors, descendants, ancestors[x], descendants[x])
        df_scores.iloc[x] = sub_score

    return df_scores
