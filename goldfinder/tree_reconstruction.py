import numpy as np
import skbio
from sklearn.metrics.pairwise import pairwise_distances
import dendropy
from dendropy.interop import raxml
from skbio import TreeNode
from six import StringIO
from ete3 import Tree
import output
import checks


class Goldfinder_Node:  # Not yet implemented, but see issue #8

    def __init__(self):
        self.index_parent = 0
        self.index_child1 = 0
        self.index_child2 = 0
        self.gene_state = 0


def main(ptree, ptree_inference, df, poutput):
    """
    ptree: path to tree newick file
    ptree_inference: method for tree inference (nj or ml)
    df: preprocessed pandas dataframe
    poutput: output path given by command line arguments
    """
    if ptree is None:
        if ptree_inference == 'nj':
            print("Constructing distance matrix")
            dm = get_dist_matrix(df)
            print("Constructing tree")
            nwk = neighbor_joining(dm)
            output.output_tree(nwk, poutput)
        else:
            print("Constructing tree")
            nwk = tree_recon_raxml(df)
            output.output_tree(nwk, poutput)
    else:
        print("Reading in phylogenetic tree")
        with open(ptree) as file:
            nwk = file.read()
        checks.tree_check(nwk, df)

    tree_struc, tip_lvlorder, branch_distances = construct_array_tree(nwk)

    return tree_struc, tip_lvlorder, branch_distances


def get_dist_matrix(df):
    """Determining the distance matrix using Jukes-Cantor correction
    df: pandas dataframe containing gene presence absence after pre-processing
    return: distance matrix
    """
    data = -0.75 * np.log(1 - pairwise_distances(df.T, metric="hamming")*(4/3))
    data[data == np.inf] = 0
    dm = skbio.DistanceMatrix(data, df.columns)

    return dm


def neighbor_joining(dist_matrix):
    """Determining the phylogenetic tree using neighbor joining
    dist_matrix: distance matrix created by function get_dist_matrix
    return: pyhlogenetic tree in newick string format
    """
    newick_str = skbio.nj(dist_matrix, disallow_negative_branch_length=False,
                          result_constructor=str)
    # setting negative branch lengths to 0
    t = Tree(newick_str)
    for node in t.traverse():
        if node.dist < 0:
            node.dist = 0

    newick_str = t.write()
    tree = TreeNode.read(StringIO(newick_str))
    tree = tree.root_at_midpoint()
    newick_str = str(tree).replace("root", "")

    return newick_str


def tree_recon_raxml(df):
    """Infering a phylogenetic tree with Maximum Likelihood using RAxML
    df: pandas presence absence table
    return: pyhlogenetic tree in newick string format
    """

    df = df.applymap(str)  # transforming into string dataframe

    seq_dict = {}
    for (columnName, columnData) in df.iteritems():
        seq_dict[columnName] = ''.join(columnData.tolist())

    pres_abs = dendropy.CharacterMatrix.from_dict(seq_dict)  # 0/1 matrix

    rx = raxml.RaxmlRunner()
    # https://github.com/amkozlov/raxml-ng/wiki/Input-data#starting-trees
    tree = rx.estimate_tree(char_matrix=pres_abs, raxml_args=['-m', 'BINCATI', '-T', '10'])
    newick_str = tree.as_string(schema="newick")

    tree = TreeNode.read(StringIO(newick_str))
    tree = str(tree.root_at_midpoint())
    newick_str = tree.replace("root", "")
    return newick_str


def construct_array_tree(nwk):
    """Transforming newick tree string into list of lists tree structure
    nwk: newick tree string
    return: list of list data structure, leaf name dictionary, list of branch lengths
    """

    # ete3 tree object for determining level order of nodes

    # check format, difference between formats
    tree = Tree(nwk, format=1)
    lvlnum = 0  # variable for labeling tree according to level order traversal

    dist = np.array([])  # storing branch length in np array for probability calculation

    # modifying tree
    for node in tree.traverse():
        if not hasattr(node, "level_order_num"):
            node.add_feature("level_order_num", lvlnum)
        # retreive branch distances and storing them in p
        dist = np.append(dist, node.dist)
        lvlnum += 1

    # index parent, index child 1, index child 2, gene presence/absence status;
    new_tree = [[0, 0, 0, 0] for i in range(lvlnum)]

    for node in tree.traverse():
        lvl = node.level_order_num

        children = node.get_children()
        count = 0
        for child in children:
            if count == 0:
                new_tree[lvl][1] = child.level_order_num
            elif count == 1:
                new_tree[lvl][2] = child.level_order_num
            else:
                break
            count += 1

        ancestors = node.get_ancestors()
        for first_anc in ancestors:
            new_tree[lvl][0] = first_anc.level_order_num
            break

    leaves = {}
    for leaf in tree.get_leaves():
        leaves[leaf.name] = leaf.level_order_num

    return new_tree, leaves, dist


def get_children_indices(tree_struc, node_index):
    """Finding indices of all children of a node
    tree_struc: tree list of lists data structure
    node_index: current index of node
    return: indices of children of node
    """
    c_ind = []

    if tree_struc[node_index][1] == 0:
        pass
        # c_ind.append(node_index)
    else:
        c_ind.append(tree_struc[node_index][1])
        c_ind.extend(get_children_indices(tree_struc, tree_struc[node_index][1]))

        c_ind.append(tree_struc[node_index][2])
        c_ind.extend(get_children_indices(tree_struc, tree_struc[node_index][2]))
    return c_ind
