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
