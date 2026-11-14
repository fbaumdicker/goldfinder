import numpy as np


def preproc(df, ppreprocess, freq_min=0.05, freq_max=0.95):
    """
    df: pandas dataframe to preprocess
    ppreprocess: boolean for removing genes with frequency outside [freq_min, freq_max]
    freq_min: minimum fraction of genomes a gene must be present in
    freq_max: maximum fraction of genomes a gene may be present in
    return: pandas dataframe
    """
    df = general_preproc(df)
    if ppreprocess:
        df = rm_insufficient_genes(df, freq_min, freq_max)
    return df


def general_preproc(df):
    """Pre-processing to remove genes present in every genome
    df: pandas dataframe containing gene presence absence from input
    return: pandas dataframe with adjusted columns
    """
    print("Preprocessing: Removing genes appearing in all genomes")
    ns = len(df.columns)  # get number of samples
    df = df.loc[(df.sum(axis=1) != ns) & (df.sum(axis=1) != 0)]  # remove genes present in every sample

    return df


def rm_insufficient_genes(df, freq_min=0.05, freq_max=0.95):
    """Pre-processing to remove genes with a frequency outside [freq_min, freq_max]
    df: pandas dataframe containing gene presence absence from input
    freq_min: minimum fraction of genomes a gene must be present in
    freq_max: maximum fraction of genomes a gene may be present in
    return: pandas dataframe with adjusted columns
    """
    print(f"Preprocessing: Removing genes with frequency below {freq_min} or above {freq_max}")
    ns = len(df.columns)  # get number of samples
    freq = df.sum(axis=1) / ns  # fraction of genomes each gene is present in
    df = df.loc[(freq >= freq_min) & (freq <= freq_max)]

    return df


def keep_gene(num_present, ns, ppreprocess, freq_min=0.05, freq_max=0.95):
    """Same filters as preproc(), applied to numbers of genomes containing a gene
    num_present: np.array (or int) of number of genomes in which a gene is present
    ns: number of genomes
    return: bool np.array, True if the gene is kept by preproc()
    """
    num_present = np.asarray(num_present)
    keep = (num_present != ns) & (num_present != 0)
    if ppreprocess:
        freq = num_present / ns
        keep = keep & (freq >= freq_min) & (freq <= freq_max)
    return keep


def filter_simulated(anc_states, desc_states, leaves_state, counts, input_leaf_mask, keep_fn):
    """Removing simulated genes that would have been removed from the input data by preproc()
    anc_states, desc_states, leaves_state, counts: simulated patterns as returned by simulation.simu
    input_leaf_mask: bool list over leaves (ordered by node index), True for tips in the input
    keep_fn: function(number of input genomes with the gene) -> bool np.array
    return: filtered anc_states, desc_states, leaves_state, counts, number of removed genes
    """
    if len(leaves_state) == 0:
        return anc_states, desc_states, leaves_state, counts, 0
    leaves = np.asarray(leaves_state)
    num_present = leaves[:, np.asarray(input_leaf_mask, dtype=bool)].sum(axis=1)
    keep = keep_fn(num_present)
    num_removed = int(np.asarray(counts)[~keep].sum())
    kept = np.nonzero(keep)[0]
    return ([anc_states[i] for i in kept], [desc_states[i] for i in kept],
            [leaves_state[i] for i in kept], [counts[i] for i in kept], num_removed)
