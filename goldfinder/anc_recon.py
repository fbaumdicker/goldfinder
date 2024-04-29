from tqdm import tqdm
import random


def recon(df, new_tree, tip_lvlorder):
    """Reconstructing the acenstral gene presence/absence state of the original input data
    using Fitch algorithm
    df: pandas dataframe containing gene presence absence from input
    new_tree: tree structure of phylogenetic tree
    tip_lvlorder: array of branch lengths
    return leaves_state: list of lists with presence/absence pattern at leaves
    return desc_state: list of lists with presence/absence pattern of children of anc
    return anc_state:
    return fitch_score_root_one: Dict {fitch_score: num trees with this score}, root state 1
    return fitch_score_root_zero: Dict {fitch_score: num trees with this score}, root state 0
    return fitch_score_root_both: Dict {fitch_score: num trees with this score}, root state 0/1
    """

    # storing gene presence/absence state for later scoring
    leaves_state = []
    desc_state = []  # representing current node
    anc_state = []  # representing ancestor of current node

    # Dictionaries {fitch_score: number of trees (genes) with this score}
    # split depending on root state
    fitch_score_root_zero = {}
    fitch_score_root_one = {}
    fitch_score_root_both = {}

    for index, row in tqdm(df.iterrows(), total=df.shape[0]):

        leaf_state = []
        descendant_state = []
        ancestor_state = []

        current_fitch_score = 0
        fitch_tree = new_tree

        # storing presence/absence of tips from data input
        row = row.to_dict()
        for key in row:
            fitch_tree[tip_lvlorder[key]][3] = row.get(key)

        # going through tree from tips to root
        for i in reversed(range(len(fitch_tree))):
            # checking if current node has a child
            if fitch_tree[i][1] != 0:
                # if same state, set state to one of the children
                if fitch_tree[fitch_tree[i][1]][3] == fitch_tree[fitch_tree[i][2]][3]:
                    fitch_tree[i][3] = fitch_tree[fitch_tree[i][1]][3]
                # if not the same, but also not "both", set state to "both"
                elif fitch_tree[fitch_tree[i][1]][3] != 2 and fitch_tree[fitch_tree[i][2]][3] != 2:
                    fitch_tree[i][3] = 2
                    current_fitch_score += 1
                # if first child is both, then use value of other child
                elif fitch_tree[fitch_tree[i][1]][3] == 2:
                    fitch_tree[i][3] = fitch_tree[fitch_tree[i][2]][3]
                # if second child is both, then use value of first child
                elif fitch_tree[fitch_tree[i][2]][3] == 2:
                    fitch_tree[i][3] = fitch_tree[fitch_tree[i][1]][3]

        # going through tree from root to tips
        for j in range(len(fitch_tree)):
            if j > 0:
                if fitch_tree[j][3] == 2:
                    fitch_tree[j][3] = fitch_tree[fitch_tree[j][1]][3]

                descendant_state.append(fitch_tree[j][3])  # store current node state
                # store node state of ancestor of current node
                ancestor_state.append(fitch_tree[fitch_tree[j][0]][3])
                # store leaf state
                if fitch_tree[j][1] == 0 and fitch_tree[j][2] == 0:
                    leaf_state.append(fitch_tree[j][3])
            # position at the root
            else:
                if fitch_tree[j][3] == 0:
                    fitch_score_root_zero[current_fitch_score] = fitch_score_root_zero.get(
                        current_fitch_score, 0) + 1

                elif fitch_tree[j][3] == 1:
                    fitch_score_root_one[current_fitch_score] = fitch_score_root_one.get(
                        current_fitch_score, 0) + 1

                elif fitch_tree[j][3] == 2:
                    fitch_tree[j][3] = random.choice([0, 1])
                    fitch_score_root_both[current_fitch_score] = fitch_score_root_both.get(
                        current_fitch_score, 0) + 1

        leaves_state.append(leaf_state)
        anc_state.append(ancestor_state)
        desc_state.append(descendant_state)
        # situation when both are both might be revisited

    return (leaves_state, desc_state, anc_state, fitch_score_root_one, fitch_score_root_zero,
            fitch_score_root_both)
