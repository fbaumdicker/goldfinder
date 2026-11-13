from tqdm import tqdm
import random


def fitch_bottom_up(tree, states):
    """First (tips to root) pass of the Fitch algorithm.
    tree: tree structure of phylogenetic tree (list of [parent, child1, child2, state]); nodes are
          in level order, i.e. every parent has a smaller index than its children
    states: list with one entry per node; tip entries must be set to 0/1. Entries of internal
            nodes are overwritten in place with their Fitch set: 0, 1 or 2 (2 = ambiguous {0, 1})
    return: Fitch score (minimum number of state changes on the tree)
    """
    fitch_score = 0
    for i in reversed(range(len(tree))):
        # checking if current node has a child
        if tree[i][1] != 0:
            state_child1 = states[tree[i][1]]
            state_child2 = states[tree[i][2]]
            # if same state, set state to one of the children
            if state_child1 == state_child2:
                states[i] = state_child1
            # if not the same, but also not "both", set state to "both"
            elif state_child1 != 2 and state_child2 != 2:
                states[i] = 2
                fitch_score += 1
            # if first child is both, then use value of other child
            elif state_child1 == 2:
                states[i] = state_child2
            # if second child is both, then use value of first child
            elif state_child2 == 2:
                states[i] = state_child1
    return fitch_score


def fitch_top_down(tree, states, root_state):
    """Second (root to tips) pass of the Fitch algorithm, assigning a concrete state to each node.
    tree: tree structure of phylogenetic tree (level order, see fitch_bottom_up)
    states: Fitch sets after fitch_bottom_up (not modified)
    root_state: concrete state (0 or 1) assigned to the root
    return leaf_state: list of states at the leaves
    return descendant_state: list of states of every non-root node
    return ancestor_state: list of states of the parent of every non-root node
    """
    final_states = list(states)
    final_states[0] = root_state

    leaf_state = []
    descendant_state = []
    ancestor_state = []

    for j in range(1, len(tree)):
        # an ambiguous node inherits the state of its parent (no change on this branch)
        if final_states[j] == 2:
            final_states[j] = final_states[tree[j][0]]

        descendant_state.append(final_states[j])  # store current node state
        # store node state of ancestor of current node
        ancestor_state.append(final_states[tree[j][0]])
        # store leaf state
        if tree[j][1] == 0 and tree[j][2] == 0:
            leaf_state.append(final_states[j])

    return leaf_state, descendant_state, ancestor_state


def recon(df, new_tree, tip_lvlorder):
    """Reconstructing the acenstral gene presence/absence state of the original input data
    using Fitch algorithm
    df: pandas dataframe containing gene presence absence from input
    new_tree: tree structure of phylogenetic tree
    tip_lvlorder: dict {tip name: node index in new_tree}
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

        # storing presence/absence of tips from data input
        states = [0] * len(new_tree)
        row = row.to_dict()
        for key in row:
            states[tip_lvlorder[key]] = row.get(key)

        current_fitch_score = fitch_bottom_up(new_tree, states)

        root_state = states[0]
        if root_state == 0:
            fitch_score_root_zero[current_fitch_score] = fitch_score_root_zero.get(
                current_fitch_score, 0) + 1

        elif root_state == 1:
            fitch_score_root_one[current_fitch_score] = fitch_score_root_one.get(
                current_fitch_score, 0) + 1

        elif root_state == 2:
            root_state = random.choice([0, 1])
            fitch_score_root_both[current_fitch_score] = fitch_score_root_both.get(
                current_fitch_score, 0) + 1

        leaf_state, descendant_state, ancestor_state = fitch_top_down(new_tree, states,
                                                                      root_state)

        leaves_state.append(leaf_state)
        anc_state.append(ancestor_state)
        desc_state.append(descendant_state)
        # situation when both are both might be revisited

    return (leaves_state, desc_state, anc_state, fitch_score_root_one, fitch_score_root_zero,
            fitch_score_root_both)


def recon_simulated(leaves_state, counts, new_tree, rng):
    """Treating simulated genes like input genes: only their tip states are kept and their
    ancestral states are reconstructed with the same Fitch algorithm (and the same rules for
    resolving ambiguity) as used for the input data in recon(). Hence the null distribution is
    computed from reconstructed histories, just like the scores of the input data.

    leaves_state: list of simulated tip states (one list per simulated pattern), with tips ordered
                  by node index as returned by simulation.simu_step
    counts: list, how often each simulated pattern occurred
    new_tree: tree structure of phylogenetic tree
    rng: random.Random instance used to choose the root state if it is ambiguous. The choice is
         made independently for every simulated gene (not once per pattern), as in recon()
    return anc_state, desc_state, leaves_state, counts: like simulation.simu, but describing the
           reconstructed histories. Genes with the same reconstructed history are merged.
    return fitch_score_root_one, fitch_score_root_zero, fitch_score_root_both: Fitch score
           distributions of the simulated genes (weighted by counts), as in recon()
    """
    leaf_nodes = [j for j in range(1, len(new_tree))
                  if new_tree[j][1] == 0 and new_tree[j][2] == 0]

    # The reconstruction only depends on the tip states, so different simulated histories
    # leading to the same tip pattern are pooled first (keeps order of first occurrence).
    tip_pattern_counts = {}
    for leaf_state, count in zip(leaves_state, counts):
        if len(leaf_state) != len(leaf_nodes):
            raise ValueError("Simulated leaf states do not match the number of tree tips")
        key = tuple(leaf_state)
        tip_pattern_counts[key] = tip_pattern_counts.get(key, 0) + count

    fitch_score_root_zero = {}
    fitch_score_root_one = {}
    fitch_score_root_both = {}

    new_anc_state = []
    new_desc_state = []
    new_leaves_state = []
    new_counts = []

    for tip_pattern, count in tqdm(tip_pattern_counts.items()):
        states = [0] * len(new_tree)
        for node, state in zip(leaf_nodes, tip_pattern):
            states[node] = state

        current_fitch_score = fitch_bottom_up(new_tree, states)

        root_state = states[0]
        if root_state == 0:
            fitch_score_root_zero[current_fitch_score] = fitch_score_root_zero.get(
                current_fitch_score, 0) + count
            root_state_counts = [(0, count)]

        elif root_state == 1:
            fitch_score_root_one[current_fitch_score] = fitch_score_root_one.get(
                current_fitch_score, 0) + count
            root_state_counts = [(1, count)]

        elif root_state == 2:
            fitch_score_root_both[current_fitch_score] = fitch_score_root_both.get(
                current_fitch_score, 0) + count
            # random root state for each of the simulated genes with this tip pattern
            num_root_zero = sum(1 for _ in range(count) if rng.choice([0, 1]) == 0)
            root_state_counts = [(0, num_root_zero), (1, count - num_root_zero)]

        for root, root_count in root_state_counts:
            if root_count == 0:
                continue
            leaf_state, descendant_state, ancestor_state = fitch_top_down(new_tree, states, root)
            new_leaves_state.append(leaf_state)
            new_desc_state.append(descendant_state)
            new_anc_state.append(ancestor_state)
            new_counts.append(root_count)

    return (new_anc_state, new_desc_state, new_leaves_state, new_counts, fitch_score_root_one,
            fitch_score_root_zero, fitch_score_root_both)
