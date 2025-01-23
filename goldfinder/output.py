import skbio
from io import StringIO
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd
import os
import clustering


def result_procedure(p_values_adj, p_values_unadj, significant_score_indices, cluster_dict,
                     clusters, locus_dict, poutput, pscore, mode, pfile_type, perform_clustering,
                     metadata, known_assoc, write_cytoscape):

    if clusters:
        print("Writing association clusters")
        cluster_file = f'{poutput}/{mode}_clusters.txt'
        write_clusters(clusters, list(p_values_adj), cluster_file, pfile_type)

        print("Preparing cluster size graphic")
        hist_file = f'{poutput}/{mode}_cluster_sizes.png'
        cluster_size_viz(clusters, hist_file)

    print("Writing significant gene pairs to output")
    gene_pair_file = f'{poutput}/{pscore}_{mode}_significant_pairs.csv'
    write_significant_gp(p_values_adj, p_values_unadj, significant_score_indices, cluster_dict,
                         locus_dict, gene_pair_file, pfile_type, perform_clustering, metadata,
                         known_assoc)

    print("Sorting output according to p-value")
    df = pd.read_csv(gene_pair_file, low_memory=False)
    sort_output(df, gene_pair_file)

    cytoscape_file = f'{poutput}/cytoscape_input.csv'
    if write_cytoscape and mode == 'association':
        print("Writing Associated Gene Pairs for Cytoscape Visualization")
        assoc_genes_cytoscape(df, cytoscape_file)

        # Calculate fraction of associated genes between clusters
        assoc_freq = clustering.dissociation_freq(cluster_dict, p_values_adj,
                                                  pfile_type in ["matrix", "tab"])

        # Also write fraction of associated genes between clusters to the cytoscape file
        clusters_assoc_cytoscape(assoc_freq, cluster_dict, poutput)

    elif write_cytoscape and mode == 'dissociation':
        print("Writing Dissociated Gene Pairs for Cytoscape Visualization")
        dissoc_genes_cytoscape(df, cytoscape_file)


def create_output_folder(poutput, pforce_output):
    """Create (and overwrite) a new output folder at location poutput
    poutput: string of path to where output should be stored
    """
    if not os.path.exists(poutput):
        os.makedirs(poutput)
    elif len(os.listdir(poutput)) and not pforce_output:
        exit('Goldfinder Error: Output directory already exists and is not empty.')


def output_tree(nwk, poutput):
    """Writing inferred phylogenetic tree as newick string and ascii representation to two txt files
    nwk: phylogenetic tree as newick string
    poutput: string of path to where output should be stored
    """
    with open(poutput + "/inf_tree_nwk.txt", "w") as text_file:
        text_file.write(nwk)
    tree = skbio.TreeNode.read(StringIO(nwk))
    with open(poutput + "/inf_tree_ascii.txt", "w") as f:
        f.write(tree.ascii_art())


def write_significant_gp(df, df_unadj, sig_indices, clusters, locus, file_name, file_type,
                         perform_clustering, metadata, known_assoc):
    """Writing significant gene pairs and their respective p-value to an output file
    df: pandas dataframe containing gene pairs and their p-value (only significant)
    sig_indices: numpy array containing indices where significant gene pairs where found
    clusters: dictionary containing genes as keys and their respective cluster as value
    file_name: name of result file, with path to output dir
    file_type: string representing input file type e.g. panx, tab or roary
    perform_clustering: bool whether clustering was performed
    metadata: user-provided metadata
    known_assoc: dict (gene1, gene2): (unadj, adj) p-val, gene1 < gene2, pair is surely associated
    """

    columns = list(df.columns)
    rows = list(df.index.values)

    known_assoc_to_write = None
    if known_assoc is not None:
        known_assoc_to_write = set(known_assoc.keys())

    with open(file_name, 'w') as f:

        # construct header based on input format
        if file_type in ["roary", "panaroo"]:
            header = ("Gene_1," + "Gene_name_1," + "Annotation_1," + "Gene_2," + "Gene_name_2," +
                      "Annotation_2," + "p-value unadj," + "p-value adj")
        elif file_type == "panx":
            header = ("Gene_1," + "Gene_name_1," + "Annotation_1," + "Locus_tags_1," + "Gene_2," +
                      "Gene_name_2," + "Annotation_2," + "Locus_tags_2," + "p-value unadj," +
                      "p-value adj")
        elif file_type in ["matrix", "tab"]:
            header = "Gene_1," + "Gene_2," + "p-value unadj," + "p-value adj"

        # add to header depending on arguments, independent of input format
        if perform_clustering:
            header += ","
            header += "Cluster"

        if metadata is not None:
            header += ","
            header += ",".join([str(col) + '_1' for col in metadata.columns])
            header += ","
            header += ",".join([str(col) + '_2' for col in metadata.columns])

        if known_assoc is not None:
            header += ","
            header += "provided_as_associated"

        header += "\n"
        f.write(header)

        # for each significant gene pair, write a line
        for x in tqdm(range(len(sig_indices[0]))):
            if file_type in ["roary", "panaroo"]:
                gene_1 = str(columns[sig_indices[0][x]]).split("/")
                gene_2 = str(rows[sig_indices[1][x]]).split("/")

                p_unadj = df_unadj.loc[str(columns[sig_indices[0][x]]),
                                       str(rows[sig_indices[1][x]])]
                p_adj = df.iloc[sig_indices[0][x], sig_indices[1][x]]

            elif file_type == "panx":
                gene_1 = str(columns[sig_indices[0][x]]).split("/")
                gene_2 = str(rows[sig_indices[1][x]]).split("/")

                p_unadj = df_unadj.loc[str(columns[sig_indices[0][x]]),
                                       str(rows[sig_indices[1][x]])]
                p_adj = df.iloc[sig_indices[0][x], sig_indices[1][x]]

            elif file_type in ["matrix", "tab"]:
                gene_1 = str(columns[sig_indices[0][x]])
                gene_2 = str(rows[sig_indices[1][x]])

                p_unadj = df_unadj.loc[str(gene_1), str(gene_2)]
                p_adj = df.iloc[sig_indices[0][x], sig_indices[1][x]]

            s, known_assoc_to_write = assemble_gp_line(gene_1, gene_2, file_type, p_unadj, p_adj,
                                                       locus, perform_clustering, clusters,
                                                       metadata, known_assoc_to_write)
            f.write(s)

        # for each non-significant but surely associated gene pair, write a line
        if known_assoc_to_write:
            # Make a shallow copy because assemble_gp_line will change the original set
            it_known_assoc = known_assoc_to_write.copy()
            for (gene_1, gene_2) in it_known_assoc:
                p_unadj, p_adj = known_assoc[(gene_1, gene_2)]
                if file_type in ["roary", "panaroo", "panX"]:
                    gene_1 = gene_1.split('/')
                    gene_2 = gene_2.split('/')
                s, known_assoc_to_write = assemble_gp_line(gene_1, gene_2, file_type, p_unadj,
                                                           p_adj, locus, perform_clustering,
                                                           clusters, metadata, known_assoc_to_write)
                f.write(s)


def assemble_gp_line(gene_1, gene_2, file_type, p_unadj, p_adj, locus_dict, perform_clustering,
                     clusters, metadata, known_assoc_to_write):
    """
    Assemble one line of the output gene pairs file

    Parameters
    ----------
    gene_1 : str
        Identifier of gene1. In case format is not tab, will contain multiple field separated by /
    gene_2 : str
        Identifier of gene2. In case format is not tab, will contain multiple field separated by /
    file_type : str
        Type of input gene absence presence matrix
    p_unadj : DataFrame
        Unadjusted p-values. Index and columns are genes. Contains NaN.
    p_adj : DataFrame
        Adjusted p-values. Index and columns are genes. Contains NaN.
    locus_dict : Dict str:str
        Present in case input is panX. Maps gene name to respective locus tag
    perform_clustering : bool
        Whether clustering was performed built
    clusters : Dict str:int
        Dict that maps gene name to its cluster number
    metadata : DataFrame
        Metadata about the genes. Gene names in index.
    known_assoc : set
        set of all known associations that still need to be written

    Returns
    -------
    s : str
        Line to write in result file
    known_assoc : set
        Updated set of all known associations that still need to be written

    """
    # this part depends on input format
    if file_type in ["roary", "panaroo"]:

        s = (gene_1[0] + "," + gene_1[1] + "," + gene_1[2] + "," + gene_2[0] + "," +
             gene_2[1] + "," + gene_2[2] + "," + str(p_unadj) + "," + str(p_adj))

        # This is the id used in the following
        gene_1 = gene_1[0]
        gene_2 = gene_2[0]

    elif file_type == "panx":

        s = (gene_1[0] + "," + gene_1[1] + "," + gene_1[2] + "," + locus_dict[gene_1[0]] + "," +
             gene_2[0] + "," + gene_2[1] + "," + gene_2[2] + "," + locus_dict[gene_2[0]] + "," +
             str(p_unadj) + "," + str(p_adj))

        # This is the id used in the following
        gene_1 = gene_1[0]
        gene_2 = gene_2[0]

    elif file_type in ["matrix", "tab"]:

        s = gene_1 + "," + gene_2 + "," + str(p_unadj) + "," + str(p_adj)

    # this part does not depend on input format but on arguments
    if perform_clustering:
        # sort_output will format this column to float if it contains "None" or "" and some integers
        cluster_name = "-"
        if (clusters.get(gene_1, 0) > 0 and clusters.get(gene_2, 0) > 0 and
                clusters[gene_1] == clusters[gene_2]):
            cluster_name = str(clusters[gene_1])

        s += ","
        s += cluster_name

    if metadata is not None:
        s += ","
        s += ",".join(metadata.loc[gene_1, :].astype(str))
        s += ","
        s += ",".join(metadata.loc[gene_2, :].astype(str))

    if known_assoc_to_write is not None:
        s += ","
        tup = (gene_1, gene_2) if gene_1 > gene_2 else (gene_2, gene_1)

        if tup in known_assoc_to_write:
            # remove the element when it is written
            known_assoc_to_write.remove(tup)
            s += 'yes'
        else:
            s += 'no'

    s += "\n"
    return s, known_assoc_to_write


def write_distribution(score_dict, poutput):

    temp_scores = list(score_dict.keys())
    temp_scores = sorted(temp_scores)

    with open(poutput, 'w') as f:
        s = "Score,Frequency" + "\n"
        f.write(s)
        for score in temp_scores:
            if score_dict[score] != 0:
                s = str(score) + "," + str(score_dict[score]) + "\n"
                f.write(s)
        f.close()


def write_log(poutput, message):
    with open(poutput+"/log.txt", 'a') as f:
        f.write(message+"\n")
        f.close()


def write_clusters(clusters, gene_names, file_name, file_type):
    """Writing for each cluster the associated genes
    clusters: mcl output consisting of clusters and the genes they contain
    gene_names: list of gene names
    """
    with open(file_name, 'w') as f:
        gene_cluster_nr = 0
        for x in tqdm(range(len(clusters))):
            if len(clusters[x]) > 1:
                cluster_size = len(clusters[x])
                gene_cluster_nr += 1
                s = ">" + str(gene_cluster_nr) + "," + str(cluster_size) + "\n"
                f.write(s)
                for gene_loc_i in range(cluster_size):
                    if file_type in ["matrix", "tab"]:
                        s = gene_names[clusters[x][gene_loc_i]]
                    else:
                        s = gene_names[clusters[x][gene_loc_i]].split("/")[0]
                    if gene_loc_i < cluster_size - 1:  # not the last gene
                        s += ','
                    f.write(s + '\n')


def cluster_size_viz(clusters, file_name):
    """Preparing a cluster size / frequency histogram
    clusters: list of lists of clusters and genes present in each cluster
    output: string representing path to output file
    """
    cluster_size = []

    for cluster in clusters:
        if len(cluster) > 1:
            cluster_size.append(len(cluster))
    cluster_size.sort()

    plt.hist(cluster_size, bins=10)
    plt.xlabel("Cluster size")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(file_name)


def cluster_dissoc(dissoc_freq, global_freq, poutput):
    """
    Writes the average dissociation of MCL clusters to file

    Parameters
    ----------
    dissoc_freq : Dict (int, int) -> float
        keys: tuple of cluster nr, values: average dissociation p-value
    global_freq : float
        Percentage of gene pairs that are significantly dissociated
    poutput : str
        User provided output folder

    Returns
    -------
    None.

    """

    print("Writing Table of Dissociation between clusters")
    with open(f'{poutput}/Dissociation_between_clusters.csv', 'w') as file:
        file.write(f'# global fraction of significant dissociations: {global_freq}\n')
        file.write('# Only contains pairs of clusters with fraction > 0\n')
        file.write('Cluster1,Cluster2,Portion of Significant Gene Pairs between Clusters\n')
        for (i, j) in dissoc_freq:
            if dissoc_freq[(i, j)] > 0:
                file.write(f'{i},{j},{dissoc_freq[(i,j)]}\n')


def assoc_genes_cytoscape(df, file_name):
    """
    Also starts the file!
    Write the internal cytoscape file used for the node layout.
    Force of associated genes is adjusted p-value - 0.5 to push strongly associated genes apart.

    Parameters
    ----------
    df : pandas DataFrame
        dataframe of associating gene pairs as constructed in this class
    file_name : str
        file name of cytoscape input

    Returns
    -------
    None.
    """

    df.rename({'Gene_1': 'Node1', 'Gene_2': 'Node2', 'p-value adj': 'Force'}, axis=1, inplace=True)
    df.drop(df.columns.difference(['Node1', 'Node2', 'Force']), axis=1, inplace=True)
    df['Force'] = df['Force'] - 0.5
    df['pair_type'] = 'gene-gene-assoc'
    df.to_csv(file_name, index=False)


def dissoc_genes_cytoscape(df, file_name):
    """
    Write to the internal cytoscape file used for the node layout.
    Force of dissociated genes is (0.5 - adjusted p-value) to pull strongly dissociated genes
    slightly together.

    Parameters
    ----------
    df : pandas DataFrame
        dataframe of dissociating gene pairs as constructed in this class
    file_name : str
        file name of cytoscape input

    Returns
    -------
    None.
    """

    df.rename({'Gene_1': 'Node1', 'Gene_2': 'Node2', 'p-value adj': 'Force'}, axis=1, inplace=True)
    df.drop(df.columns.difference(['Node1', 'Node2', 'Force']), axis=1, inplace=True)
    df['Force'] = 0.5 - df['Force']
    df['pair_type'] = 'gene-gene-dissoc'
    df.to_csv(file_name, mode='a', header=False, index=False)


def clusters_cytoscape(dissoc_freq, cluster_dict, poutput):
    """
    Write to the internal cytoscape file used for the node layout.
    Force of Genes to their cluster nodes is 5 so these are grouped together for sure.
    Force between clusters is the fraction of dissociated genes between two clusters + 1 to place
    dissociated clusters close together. This way, clusters with any connection (assoc or dissoc)
    will be placed closeby.

    Parameters
    ----------
    dissoc_freq : Dict (int, int) -> float
        keys: tuple of cluster nr, values: average dissociation p-value
    cluster_dict : Dict str -> int
        Dictionary of gene_name to ID of MCL cluster, which was generated using gene associations
    poutput : str
        User provided output folder

    Returns
    -------
    None.
    """

    print("Writing Cluster Nodes for Cytoscape Visualization")
    with open(f'{poutput}/cytoscape_input.csv', 'a') as file:

        # write gene-cluster pairs
        for gene in cluster_dict:
            file.write(f'{gene},cl_{cluster_dict[gene]},5,gene-cluster-member\n')

        # write cluster-cluster-dissoc pairs
        for (i, j) in dissoc_freq:
            if dissoc_freq[(i, j)] > 0:
                file.write(f'cl_{i},cl_{j},{dissoc_freq[(i,j)] + 1},cluster-cluster-dissoc\n')


def clusters_assoc_cytoscape(assoc_freq, cluster_dict, poutput):
    """
    Write to the internal cytoscape file used for the node layout.
    Force between clusters is the fraction of associated genes between two clusters + 1 to place
    associated clusters closer together. This way, clusters with any connection (assoc or dissoc)
    will be placed closeby.

    Parameters
    ----------
    assoc_freq : Dict (int, int) -> float
        keys: tuple of cluster nr, values: average association p-value
    cluster_dict : Dict str -> int
        Dictionary of gene_name to ID of MCL cluster, which was generated using gene associations
    poutput : str
        User provided output folder

    Returns
    -------
    None.
    """

    print("Writing Cluster Nodes for Cytoscape Visualization")
    with open(f'{poutput}/cytoscape_input.csv', 'a') as file:

        # write cluster-cluster-assoc pairs
        for (i, j) in assoc_freq:
            if assoc_freq[(i, j)] > 0:
                file.write(f'cl_{i},cl_{j},{assoc_freq[(i,j)] + 1},cluster-cluster-assoc\n')


def write_adjac_matrix(adj_mtx, poutput):
    """
    Write the adjacency matrix used for clustering to file

    Parameters
    ----------
    adj_mtx : Numpy Array
        This adjacency matrix is generated in clsutering.adjac_matrix
    poutput : str
        argument provided by the user. Output directory.

    Returns
    -------
    None.

    """
    pd.DataFrame(adj_mtx).to_csv(poutput + '/adjacency_matrix.csv')


def sort_output(df, file_name):
    """Sorting the output of significant gene pairs according to their adjusted p-value
    """

    df.sort_values(by=['p-value adj', 'Gene_1']).to_csv(file_name, index=False)
