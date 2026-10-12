def preproc(df, ppreprocess):
    """
    df: pandas dataframe to preprocess
    ppreprocess: boolean for removing genes present in less than 5% of genomes
    return: pandas dataframe
    """
    df = general_preproc(df)
    if ppreprocess:
        df = rm_insufficient_genes(df)
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


def rm_insufficient_genes(df):
    """Pre-processing to remove genes present in less than 5% of genomes
    df: pandas dataframe containing gene presence absence from input
    return: pandas dataframe with adjusted columns
    """
    print("Preprocessing: Removing genes appearing only in 5% of genomes")
    ns = len(df.columns)  # get number of samples
    tresh = round(ns * 0.05)  # int representing 5% of number of samples
    df = df.loc[df.sum(axis=1) >= tresh]  # remove genes present in less than 5% of sample

    return df
