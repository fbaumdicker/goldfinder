import markov_clustering as mc
import numpy as np
# import output
from itertools import combinations
from tqdm import tqdm


def cluster_procedure(p_values_adj, sig_lvl, pinflation, poutput, perform_clustering):
    """
    return:
        dict_cluster: dictionary of gene_name: 1-based cluster number
        clusters: list of tuples of indices of genes in the cluster represented by tuple
    """

    if perform_clustering:
        print("Constructing association networks")
        adj_mtx = adjac_matrix(p_values_adj, sig_lvl)
        # output.write_adjac_matrix(adj_mtx, poutput) # Currently not needed
        clusters = mcl(adj_mtx, pinflation)
        clusters = get_relevant_clusters(clusters)
        dict_cluster = cluster_dict(clusters, list(p_values_adj.columns))
    else:
        clusters = list()
        dict_cluster = dict()

    return dict_cluster, clusters


def dissociation_freq(cluster_dict, p_values_adj, simple_idx):
    """
    Calculate average dissociation between MCL clusters, which were calculated based on association
    between genes.
    If p_values_adj are p_values for association, this can calculate average association between
    MCL clusters just as well.

    Parameters
    ----------
    cluster_dict : Dict str -> int
        Dictionary of gene_name to ID of MCL cluster, which was generated using gene associations
    p_values_adj : DataFrame
        p-values of dissociation of genes
    simple_idx : bool
        Whether index of central dataframe is simple or contains additional information. Depends on
        input format.

    Returns
    -------
    dissoc_freq : Dict (int, int) -> float
        keys: tuple of cluster IDs. values: percentage of significantly dissociated gene pairs
        between the two clusters
    """

    # Reverse cluster_dict to cluster_ID -> [List of gene_names]
    rev_dict = {}
    for gene_name in cluster_dict:
        rev_dict[cluster_dict[gene_name]] = rev_dict.get(cluster_dict[gene_name], []) + [gene_name]

    dissoc_freq = {}

    # p_values_adj's index might contain additional information, separated by /
    id_to_index = {}
    if simple_idx:
        all_ids = p_values_adj.index
    else:
        for idx in p_values_adj.index:
            id_to_index[idx.split('/')[0]] = idx
        all_ids = set(id_to_index.keys())

    if len(rev_dict) > 1:
        # for each pair of clusters cl1, cl2
        for cl1, cl2 in tqdm(combinations(rev_dict.keys(), 2),
                             total=int(len(rev_dict) * (len(rev_dict) - 1) / 2)):

            if simple_idx:
                indices_1 = all_ids.intersection(rev_dict[cl1])
                indices_2 = all_ids.intersection(rev_dict[cl2])
            else:
                indices_1 = [id_to_index[iden] for iden in all_ids.intersection(rev_dict[cl1])]
                indices_2 = [id_to_index[iden] for iden in all_ids.intersection(rev_dict[cl2])]

            # divisor will never be 0 because clusters dict has no empty clusters
            dissoc_freq[(cl1, cl2)] = (p_values_adj.loc[indices_1, indices_2]
                                       .notna().to_numpy().sum()
                                       / (len(rev_dict[cl1]) * len(rev_dict[cl2])))
    else:
        print("Warning: less than 2 clusters present!")

    return dissoc_freq


def adjac_matrix(df, alpha):
    """Building of adjacency matrix
    df: data frame containing p values
    alpha: significance level
    return:
    """

    df = df[(df < float(alpha))]

    np.fill_diagonal(df.values, np.nan)

    df = df.dropna(axis=0, how="all")

    df = df.dropna(axis=1, how="all")

    df = df.notnull().astype('int')

    np.fill_diagonal(df.values, 0)

    return df.to_numpy()


def mcl(adj_mat, inf):
    """Mcl algorithm
    adj_mat: adjacency matrix filled with ones (connection) and zeros (no connection) between genes
    inf: float representing the inflation point (influencing resolution of clustering)
    return:
    """
    result = mc.run_mcl(adj_mat, inflation=inf)
    clusters = mc.get_clusters(result)
    return clusters


def get_relevant_clusters(clusters):

    relev_clusters = []
    for cluster in clusters:
        if len(cluster) > 1:
            relev_clusters.append(cluster)
    return relev_clusters


def cluster_dict(clusters, gene_names):
    """Representing genes and their respective clusters as dictionary
    clusters: clusters containing genes as returned by mcl
    gene_names: list containing gene names
    return: dictionary of gene_name: 1-based cluster number
    """
    cl_dict = {}

    for x in range(len(clusters)):
        for gene_loc in clusters[x]:
            # set cluster number of curr gene as 1-based
            cl_dict[gene_names[gene_loc].split("/")[0]] = x+1

    return cl_dict
