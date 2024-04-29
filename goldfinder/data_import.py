import pandas as pd
import numpy as np
import json
from Bio import SeqIO
from tqdm import tqdm


def load_input(pinput, filetype, pmetadata, pknown_associations):
    """
    pinput: string path to input roary/panx/tab-delimited file
    filetype: string roary, tab or panx
    pmetadata: input parameter metadata
    pknown_associations: input parameter known_associations
    return: df: pandas dataframe
    return: locus_dict: in case of panx input file, additional annotations are stored here
    """
    locus_dict = 0
    if filetype == "roary":
        df = load_roary(pinput)
    elif filetype == "tab":
        df = load_tab_delimited(pinput)
    elif filetype == "panx":
        df, locus_dict = load_panx(pinput)

    metadata = None
    if pmetadata is not None:
        if filetype == "tab":
            metadata = load_metadata(pmetadata, df.index)
        else:  # TODO not sure this also applies to panx case
            metadata = load_metadata(pmetadata, pd.Index([i.split('/')[0] for i in df.index]))

    num_known_assoc = 0
    known_assoc = None
    if pknown_associations is not None:
        try:  # This can either be an int
            num_known_assoc = int(pknown_associations)
            if num_known_assoc > (df.shape[0] * (df.shape[0] - 1)):
                exit('Number of known associations is larger than number of gene pairs.')
            if num_known_assoc < 0:
                exit('Number of known associations was parsed as negative.')
        except ValueError:  # or a file
            num_known_assoc, known_assoc = load_known_assoc(pknown_associations, df.index,
                                                            filetype == "tab")

    return df, locus_dict, metadata, num_known_assoc, known_assoc


def load_known_assoc(pknown_assoc, full_gene_info, filetype_tab):
    df = pd.read_csv(pknown_assoc, sep='\t', header=None)

    if df.shape[1] != 2:
        exit('List of known associations does not contain exactly 2 columns')

    # filter duplicates and same-gene-pairs
    num_pairs = df.shape[0]
    df = df.loc[(~df.duplicated()) & (df.iloc[:, 0] != df.iloc[:, 1]), :]

    if df.shape[0] != num_pairs:
        print(f'Warning: List of known associations included {num_pairs - df.shape[0]} duplicated '
              'or same-gene-pairs. These were removed.')
        num_pairs = df.shape[0]

    # Dict gene_id: full_gene_info for this function. In case of tab format, both are identical
    genes_dict = {}
    for gene_info in full_gene_info:
        genes_dict[gene_info.split('/')[0]] = gene_info

    # Check whether data contains a header, if yes remove it
    if not df.iloc[0, :].isin(genes_dict.keys()).all():
        df = df.iloc[1:, :]

    # Check whether all entries in the list are genes
    if not df.isin(genes_dict.keys()).all().all():
        err_msg = 'List of known associations contains entries that are not genes, e.g. '
        err_msg += ', '.join(list(set(df.stack()).difference(genes_dict.keys()))[:3])
        exit(err_msg)

    # dict (gene1, gene2): (Nan, Nan) with gene1 > gene2 lexicographically. Nans will become (unadj,
    # adj) p-values. In case of format not tab, keys will match gene_abs_pres dataframe index
    known_assoc = {}
    for tup in (tuple(sorted(genes, reverse=True)) for genes in df.values):
        if filetype_tab:
            known_assoc[tup] = (np.nan, np.nan)
        else:
            known_assoc[tuple(genes_dict[i] for i in tup)] = (np.nan, np.nan)

    return num_pairs, known_assoc


def load_roary(roary):
    """Load gene presence absence csv file from Roary into a pandas dataframe
    roary: roary csv file
    return: pandas dataframe
    """
    df = pd.read_csv(roary, index_col=0, low_memory=False)

    annotation = df.iloc[:, 0:2]
    df = df.drop(df.iloc[:, 0:13], axis=1)  # remaining columns are samples/species
    # df = df.notnull().astype('int') #empty cells'value replaced with 0 and non-empty with 1

    # incorporating annotation in gene_name
    tmp = {}
    for index, row in annotation.iterrows():
        tmp[index] = str(index)+"/"+str(row.iloc[0])+"/"+str(row.iloc[1]).replace(",", "")
    # renaming gene name to gene name + annotation
    df.rename(index=tmp, inplace=True)

    # convert gene name (or x in small example) to 1, and nan to 0
    df = df.map(lambda x: 1 if type(x) == str else 0)

    return df


def load_tab_delimited(tab_del):
    """Load tab-delimited gene presence absence file into a pandas dataframe
    tab_del: tab-delimited text file
    return: pandas dataframe
    """
    df = pd.read_csv(tab_del, sep="\t", header=None, low_memory=False)
    genes = df.iloc[:, 0].unique()  # row names
    species = df.iloc[:, 1].unique()  # column names

    in_matrix = pd.DataFrame(np.zeros((len(genes), len(species))), index=genes, columns=species)
    for index, row in df.iterrows():
        in_matrix.loc[row[0], row[1]] = 1

    return in_matrix.astype(int)


def load_panx(panx_folder):
    """Load panx output folder and convert it into a pandas dataframe
    panx_folder: panx output folder
    return: pandas dataframe
    """

    tab_del = []

    locus = {}
    print('Reading Gene Information')
    with open(panx_folder+"/vis/geneCluster.json") as json_file:
        data = json.load(json_file)

        for obj in tqdm(data):
            cluster_name = obj['msa']
            ann = obj['ann'].replace(",", "")
            gene_name = obj['GName']

            locus[cluster_name] = obj['locus']

            for record in SeqIO.parse(panx_folder+"/geneCluster/"+cluster_name+".faa", "fasta"):
                tab_del.append([cluster_name+"/"+gene_name+"/"+ann, str(record.id).split("|")[0]])

    tmp = np.asarray(tab_del)

    genes = np.unique(tmp[:, 0])  # row names
    species = np.unique(tmp[:, 1])  # column names

    # presence absence matrix filled with zeros
    in_matrix = pd.DataFrame(np.zeros((len(genes), len(species))), index=genes, columns=species)
    # fill matrix with presence
    print('Reading Gene Presence Abscence')
    for row in tqdm(tmp):
        in_matrix.loc[row[0], row[1]] = 1

    return in_matrix.astype(int), locus


def load_metadata(pmetadata, gene_names):

    # Read metadata
    metadata = pd.read_csv(pmetadata, sep='\t', index_col=0, header=0, dtype=str)

    # Check gene number, identification
    if len(metadata.index.intersection(gene_names)) == 0:
        exit("None of the genes in metadata table match the genes in presence-absence table.")

    if sorted(list(metadata.index)) != sorted(list(gene_names)):
        print("Warning: Genes in metadata and presence-absence table do not match exactly.")

    return metadata
