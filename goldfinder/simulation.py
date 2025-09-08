from tqdm import tqdm
import random
import numpy as np
import output
import tree_reconstruction


def simu(dist_one, dist_zero, dist_both, new_tree, branch_dist, s):
    """Simulating genes
    dist_one: homoplasy distribution for root state 1
    dist_zero: homoplasy distribution for root state 0
    dist_both: homoplasy distribution for root state 0/1
    new_tree: tree data structure
    branch_dist: list containing tree branch distances
    s: number of genes to simulate
    """
    p_norm = branch_dist/branch_dist.sum()  # normalizing branch length for probability

    # Determine how many genes are allocated to each distribution
    total = sum(dist_both.values()) + sum(dist_one.values()) + sum(dist_zero.values())
    distribution_dict = {}
    distribution_dict["1"] = round(sum(dist_one.values())/total * s)
    distribution_dict["0"] = round(sum(dist_zero.values())/total * s)
    distribution_dict["2"] = round(sum(dist_both.values())/total * s)

    # if for whatever reason they do not add up to s, reduce/increase homoplasy dist with most genes
    if sum(distribution_dict.values()) != s:
        m = max(distribution_dict, key=distribution_dict.get)
        if sum(distribution_dict.values()) < s:
            distribution_dict[str(m)] = (distribution_dict[str(m)] +
                                         (s - sum(distribution_dict.values())))
        else:
            distribution_dict[str(m)] = (distribution_dict[str(m)] -
                                         (sum(distribution_dict.values())-s))

    patterns = []
    # anc_states: list of (list of gene states of all ancestors) for each simulation
    anc_states = []
    desc_states = []
    leaves = []
    # counts: How often does a pattern occur? Same ordering as the patterns list
    counts = []

    if distribution_dict["1"] != 0:

        dist_one = [dist_one.get(i+1, 0) for i in range(max(dist_one.keys()))]
        dist_one = [round(x/sum(dist_one)*distribution_dict["1"]) for x in dist_one]

        if sum(dist_one) != distribution_dict["1"]:
            m = dist_one.index(max(dist_one))
            if sum(dist_one) < s:
                dist_one[m] = dist_one[m] + (distribution_dict["1"] - sum(dist_one))
            else:
                dist_one[m] = dist_one[m] - (sum(dist_one)-distribution_dict["1"])

        n_mts_one = []
        for j in range(len(dist_one)):
            n_mts_one.extend([j+1]*dist_one[j])
        random.shuffle(n_mts_one)

        # for each gene determine which branch shall be mutated
        branch_mutate_one = [np.random.choice(len(branch_dist), size=int(x), replace=False,
                                              p=p_norm) for x in n_mts_one]
        branch_mutate_one = [sorted(x.tolist()) for x in branch_mutate_one]

        simu_step(1, branch_mutate_one, patterns, new_tree, counts, leaves, desc_states, anc_states)

    if distribution_dict["0"] != 0:
        dist_zero = [dist_zero.get(i+1, 0) for i in range(max(dist_zero.keys()))]
        dist_zero = [round(x/sum(dist_zero)*distribution_dict["0"]) for x in dist_zero]

        if sum(dist_zero) != distribution_dict["0"]:
            m = dist_zero.index(max(dist_zero))
            if sum(dist_zero) < s:
                dist_zero[m] = dist_zero[m] + (distribution_dict["0"] - sum(dist_zero))
            else:
                dist_zero[m] = dist_zero[m] - (sum(dist_zero)-distribution_dict["0"])

        n_mts_zero = []
        for j in range(len(dist_zero)):
            n_mts_zero.extend([j+1]*dist_zero[j])
        random.shuffle(n_mts_zero)

        # for each gene determine which branch shall be mutated
        branch_mutate_zero = [np.random.choice(len(branch_dist), size=int(x), replace=False,
                                               p=p_norm) for x in n_mts_zero]
        branch_mutate_zero = [sorted(x.tolist()) for x in branch_mutate_zero]

        simu_step(0, branch_mutate_zero, patterns, new_tree,
                  counts, leaves, desc_states, anc_states)

    if distribution_dict["2"] != 0:
        dist_both = [dist_both.get(i+1, 0) for i in range(max(dist_both.keys()))]
        dist_both = [round(x/sum(dist_both)*distribution_dict["2"]) for x in dist_both]

        if sum(dist_both) != distribution_dict["2"]:
            m = dist_both.index(max(dist_both))
            if sum(dist_both) < s:
                dist_both[m] = dist_both[m] + (distribution_dict["2"] - sum(dist_both))
            else:
                dist_both[m] = dist_both[m] - (sum(dist_both)-distribution_dict["2"])

        n_mts_both = []
        for j in range(len(dist_both)):
            n_mts_both.extend([j+1]*dist_both[j])
        random.shuffle(n_mts_both)

        # for each gene determine which branch shall be mutated
        branch_mutate_both = [np.random.choice(len(branch_dist), size=int(x), replace=False,
                                               p=p_norm) for x in n_mts_both]
        branch_mutate_both = [sorted(x.tolist()) for x in branch_mutate_both]

        simu_step(2, branch_mutate_both, patterns, new_tree,
                  counts, leaves, desc_states, anc_states)

    return anc_states, desc_states, leaves, counts, dist_one, dist_zero, dist_both


def simu_step(dist_name, branch_mutate, patterns, new_tree, counts, leaves_state, desc_state,
              anc_state):
    """Performing state changes for tree nodes based on which branches were selected
    """
    root = dist_name
    if root == 2:
        root = random.choice([0, 1])

    for x in tqdm(branch_mutate):

        simu_tree = list(new_tree)

        if [root, x] not in patterns:
            leaf_state = []
            descendant_state = []
            ancestor_state = []

            # setting every node to root state
            for node in range(len(simu_tree)):
                simu_tree[node][3] = root

            # mutate nodes
            # TODO this is inefficient: for each mutation, the whole subtree's genotype is changed.
            # If a mutation occurs in the subtree, nodes' genotypes are changed multiple times
            for node in x:
                if simu_tree[simu_tree[node][0]][3] == 1:
                    simu_tree[node][3] = 0
                else:
                    simu_tree[node][3] = 1
                # modifying child nodes
                children = tree_reconstruction.get_children_indices(simu_tree, node)
                for child in children:
                    simu_tree[child][3] = simu_tree[node][3]

            # collect states in lists
            for node in range(len(simu_tree)):
                if node > 0:
                    descendant_state.append(simu_tree[node][3])  # store current node state
                    # store node state of ancestor of current node
                    ancestor_state.append(simu_tree[simu_tree[node][0]][3])
                    # store leaf state
                    if simu_tree[node][1] == 0 and simu_tree[node][2] == 0:
                        leaf_state.append(simu_tree[node][3])

            patterns.append([root, x])
            counts.append(1)
            leaves_state.append(leaf_state)
            anc_state.append(ancestor_state)
            desc_state.append(descendant_state)

        else:
            counts[patterns.index([root, x])] += 1


def additional_simu_new(scores, dist_one, dist_zero, dist_both, new_tree, branch_dist, s,
                        poutput, padditional):

    p_norm = branch_dist/branch_dist.sum()  # normalizing branch length for probability
    patterns = []
    anc_states = []
    desc_states = []
    leaves = []
    counts = []

    # transforming scores dictionary into sorted score array
    score_list = sorted([key for key in scores if scores[key] != 0])
    score_frequ = [scores[x] for x in score_list]

    # Calculate percentiles for score distribution
    cumsum = np.cumsum(score_frequ)
    percentiles = cumsum/cumsum[-1]

    # Get score of user defined percentile
    k_index = np.where(percentiles > padditional)[0][0]
    k = score_list[k_index]

    k_frac = 0
    dist_sum = 0

    # From here on, the score limit k is used with Fitch scores. This is by design!
    # Fitch score = num mutations on tree, min(Fitch(gene1), Fitch(gene2)) >= simultaneous score
    for dist in [dist_one, dist_zero, dist_both]:
        if len(dist) > 0:
            dist_sum += sum(dist.values())
            for score in dist.keys():
                if score >= k:
                    k_frac += dist[score]

    k_frac = k_frac / dist_sum

    output.write_log(poutput, f"Threshold k of the more extreme simultaneous scores: {k}")
    output.write_log(poutput, f"Fraction k_frac of simultaneous scores in extreme region: {k_frac}")

    # Stop additional simulation if it is not more fine-grained than default simulation
    if k_frac == 1:
        return None, None, None, None, k, k_frac

    # Determine how many genes are allocated to each distribution
    total = sum(dist_both.values()) + sum(dist_one.values()) + sum(dist_zero.values())
    num_sim_genes_with_root_state = {}
    num_sim_genes_with_root_state["1"] = round(sum(dist_one.values())/total * s)
    num_sim_genes_with_root_state["0"] = round(sum(dist_zero.values())/total * s)
    num_sim_genes_with_root_state["2"] = round(sum(dist_both.values())/total * s)

    # if for whatever reason they do not add up to s, reduce/increase homoplasy dist with most genes
    if sum(num_sim_genes_with_root_state.values()) != s:
        m = max(num_sim_genes_with_root_state, key=num_sim_genes_with_root_state.get)
        if sum(num_sim_genes_with_root_state.values()) < s:
            num_sim_genes_with_root_state[str(m)] = num_sim_genes_with_root_state[str(m)] + \
                (s - sum(num_sim_genes_with_root_state.values()))
        else:
            num_sim_genes_with_root_state[str(m)] = num_sim_genes_with_root_state[str(m)] - \
                (sum(num_sim_genes_with_root_state.values())-s)

    if num_sim_genes_with_root_state["1"] != 0:

        # Only look at Fitch scores >= k
        d = {key: dist_one[key] for key in dist_one if key >= k}
        sum_values = sum(d.values())

        if sum_values > 0:

            # Normalize and multiply by share of s that trees with root state 1 take up
            d = {key: round(d[key] / sum_values * num_sim_genes_with_root_state["1"])
                 for key in d}

            # If values do not add up to desired number, change max by the difference
            if sum(d.values()) != num_sim_genes_with_root_state["1"]:
                m = max(d, key=lambda x: d[x])
                d[m] = d[m] + (num_sim_genes_with_root_state["1"] - sum(d.values()))

            # n_mutations: Number of mutations for trees with root state 1
            # will contain fitch scores, and each one as often as relatively appropriate from input
            # distribution and given number of genes to be simulated
            n_mutations = []
            for key in d:
                n_mutations.extend([key] * d[key])
            random.shuffle(n_mutations)

            # for each gene determine which branch shall be mutated
            # Fitch score determines how many branches are mutated
            branch_mutate_one = [np.random.choice(len(branch_dist), size=int(x), replace=False,
                                                  p=p_norm) for x in n_mutations]
            branch_mutate_one = [sorted(x.tolist()) for x in branch_mutate_one]

            simu_step(1, branch_mutate_one, patterns, new_tree,
                      counts, leaves, desc_states, anc_states)

    if num_sim_genes_with_root_state["0"] != 0:

        # Only look at Fitch scores >= k
        d = {key: dist_zero[key] for key in dist_zero if key >= k}
        sum_values = sum(d.values())

        if sum_values > 0:

            # Normalize and multiply by share of s that trees with root state 0 take up
            d = {key: round(d[key] / sum_values * num_sim_genes_with_root_state["0"])
                 for key in d}

            # If values do not add up to desired number, change max by the difference
            if sum(d.values()) != num_sim_genes_with_root_state["0"]:
                m = max(d, key=lambda x: d[x])
                d[m] = d[m] + (num_sim_genes_with_root_state["0"] - sum(d.values()))

            n_mutations = []
            for key in d:
                n_mutations.extend([key] * d[key])
            random.shuffle(n_mutations)

            # for each gene determine which branch shall be mutated
            branch_mutate_zero = [np.random.choice(len(branch_dist), size=int(x), replace=False,
                                                   p=p_norm) for x in n_mutations]
            branch_mutate_zero = [sorted(x.tolist()) for x in branch_mutate_zero]

            simu_step(0, branch_mutate_zero, patterns, new_tree,
                      counts, leaves, desc_states, anc_states)

    if num_sim_genes_with_root_state["2"] != 0:

        # Only look at Fitch scores >= k
        d = {key: dist_both[key] for key in dist_both if key >= k}
        sum_values = sum(d.values())

        if sum_values > 0:

            # Normalize and multiply by share of s that trees with root state 1/0 take up
            d = {key: round(d[key] / sum_values * num_sim_genes_with_root_state["2"])
                 for key in d}

            # If values do not add up to desired number, change max by the difference
            if sum(d.values()) != num_sim_genes_with_root_state["2"]:
                m = max(d, key=lambda x: d[x])
                d[m] = d[m] + (num_sim_genes_with_root_state["2"] - sum(d.values()))

            n_mutations = []
            for key in d:
                n_mutations.extend([key] * d[key])
            random.shuffle(n_mutations)

            # for each gene determine which branch shall be mutated
            branch_mutate_both = [np.random.choice(len(branch_dist), size=int(x), replace=False,
                                                   p=p_norm) for x in n_mutations]
            branch_mutate_both = [sorted(x.tolist()) for x in branch_mutate_both]

            simu_step(2, branch_mutate_both, patterns, new_tree,
                      counts, leaves, desc_states, anc_states)

    return anc_states, desc_states, leaves, counts, k, k_frac
