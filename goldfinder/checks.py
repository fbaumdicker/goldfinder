import os
import pandas as pd
from ete3 import Tree


def check_input_file(pinput, pfiletype):
    """Checking input files
    pinput: path to gene presence/absence file
    filetype: string describing if roary, tab-delimited or panx
    """

    if pfiletype == "roary":
        roary_check(pinput)
    elif pfiletype == "panx":
        panx_check(pinput)


def tree_check(nwk, df):
    """Checking if tree leaves and input file genomes match
    nwk: newick string read from user provided file
    df: gene absence presence data as pandas dataframe
    """

    tree_struc = Tree(nwk, format=1)
    leaf_names = {leaf.name for leaf in tree_struc}
    if set(df.columns).issubset(leaf_names):
        print("Tree coincides with input file")
    else:
        exit("Leaf names of provided tree do not match with species name of provided input file")


def roary_check(roary_file):
    """Checking if input file adheres to roary file type
    roary_file: path to input roary file
    """
    try:
        df = pd.read_csv(roary_file, index_col=0, low_memory=False)
    except PermissionError:
        exit('Encountered permission error. Did you specify the correct input format using -f?')
    except FileNotFoundError:
        exit('Could not find specified file. Perhaps a typo or you gave a panX folder without '
             'specifying -f panx?')

    if list(df.columns.values)[0:13] != ['Non-unique Gene name', 'Annotation', 'No. isolates',
                                         'No. sequences', 'Avg sequences per isolate',
                                         'Genome Fragment', 'Order within Fragment',
                                         'Accessory Fragment', 'Accessory Order with Fragment',
                                         'QC', 'Min group size nuc', 'Max group size nuc',
                                         'Avg group size nuc']:
        print("File does not seem to be in roary format")
        exit()
    else:
        print("File is in roary format.")


def panx_check(panx_folder):
    """Checking if the input panx folder contains all of the necessary files
    panx_folder: path to input panx folder
    """
    if os.path.isdir(panx_folder):
        print("Input directory exists")

        if os.path.isfile(panx_folder+"/vis/geneCluster.json"):
            print("File is there")
        else:
            exit("geneCluster.json could not be found")

        if os.path.isdir(panx_folder+"/geneCluster"):
            print("geneCluster directory exists")
            faa_exists = False
            for fname in os.listdir(panx_folder+"/geneCluster"):
                if fname.endswith(".faa"):
                    print("There exists at least one file ending in .faa")
                    faa_exists = True
                    break
            if not faa_exists:
                exit("No file with .faa was found.")

    else:
        exit("Directory does not exist. IOError has occured.")
