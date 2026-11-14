"""Calibrating the number of mutations placed on the tree for the simulated (null) genes.

The simulation places M mutations on randomly chosen branches. Previously M was set to the Fitch
score s of an input gene. However, Fitch reconstruction can only remove changes (s <= M): hidden
reversals and parallel changes are not identifiable from the tips. Simulated genes with M = s thus
have a lower Fitch score after reconstruction than the input genes they are meant to represent.

Lookup table approach: For M = 1, ..., M_max, place M mutations with exactly the same procedure as
simulation.simu() many times and record the Fitch score s of the resulting tip pattern, giving
P(s | M). By Bayes' rule P(M | s) ~ P(s | M) * prior(M). Each simulated gene that represents an
input gene with Fitch score s draws its number of mutations from P(M | s) (the distribution, not
the expected value).

The prior over M is estimated from the input Fitch score distribution with the EM algorithm
(iterative Bayesian unfolding): starting from a uniform prior over 1, ..., M_max (which is exactly
the one-step lookup table), each iteration sets the prior to the distribution of M implied by the
current P(M | s). At convergence, simulated genes reconstructed with Fitch reproduce the input
Fitch score distribution.

If simulated genes are filtered like the input genes (constant genes and, with --preprocess,
genes outside the frequency limits), the lookup table records P(s, gene passes filter | M) and the
filtered genes are treated as unobserved, so the genes surviving the filter match the input genes.
"""
import numpy as np
from tqdm import tqdm


def tree_arrays(new_tree):
    """Parent and child indices of the tree structure as numpy arrays (level order is checked)
    new_tree: tree structure (list of [parent, child1, child2, state])
    return: parent, child1, child2 (np.arrays), leaf_nodes (list of leaf node indices)
    """
    parent = np.array([node[0] for node in new_tree])
    child1 = np.array([node[1] for node in new_tree])
    child2 = np.array([node[2] for node in new_tree])
    if np.any(parent[1:] >= np.arange(1, len(new_tree))):
        raise ValueError("Tree nodes are not in level order (parent index must be smaller)")
    leaf_nodes = [j for j in range(1, len(new_tree)) if child1[j] == 0 and child2[j] == 0]
    return parent, child1, child2, leaf_nodes


def simulate_mutation_sets(rng, branch_dist, num_mutations, replicates):
    """Placing num_mutations mutations on distinct branches, chosen with probability proportional
    to branch length, exactly as in simulation.simu()
    rng: np.random.RandomState
    return: list of np.arrays of mutated node indices
    """
    p_norm = branch_dist / branch_dist.sum()
    return [rng.choice(len(branch_dist), size=int(num_mutations), replace=False, p=p_norm)
            for _ in range(replicates)]


def tips_and_fitch(tree_arr, mutation_sets):
    """Tip states (for root state 0) and Fitch scores of genes with given mutated nodes.
    Vectorized equivalent of simulation.simu_step followed by anc_recon.fitch_bottom_up: the state
    of a node is the root state XOR the parity of mutated nodes on the path from the root to the
    node (including the root and the node itself).
    tree_arr: output of tree_arrays()
    mutation_sets: list of arrays of mutated node indices, one per gene
    return: leaf states (bool array, leaves x genes, leaves ordered by node index), Fitch scores
    """
    parent, child1, child2, leaf_nodes = tree_arr
    num_nodes = len(parent)
    num_genes = len(mutation_sets)

    mutated = np.zeros((num_nodes, num_genes), dtype=bool)
    rows = np.concatenate([np.asarray(x, dtype=np.int64) for x in mutation_sets])
    cols = np.repeat(np.arange(num_genes), [len(x) for x in mutation_sets])
    mutated[rows, cols] = True

    state = np.zeros((num_nodes, num_genes), dtype=bool)
    state[0] = mutated[0]
    for v in range(1, num_nodes):
        state[v] = state[parent[v]] ^ mutated[v]

    # Fitch sets as bits: state 0 -> 1, state 1 -> 2, ambiguous {0, 1} -> 3. Internal nodes start
    # with state 0, as in anc_recon (only matters for malformed trees with single children).
    fitch_sets = np.ones((num_nodes, num_genes), dtype=np.uint8)
    fitch_sets[leaf_nodes] = state[leaf_nodes].astype(np.uint8) + 1
    scores = np.zeros(num_genes, dtype=np.int64)
    for i in reversed(range(num_nodes)):
        if child1[i] != 0:
            set1 = fitch_sets[child1[i]]
            set2 = fitch_sets[child2[i]]
            intersection = set1 & set2
            empty = intersection == 0
            fitch_sets[i] = np.where(empty, set1 | set2, intersection)
            scores += empty
    return state[leaf_nodes], scores


def fitch_lookup_table(new_tree, branch_dist, input_leaf_mask, keep_fn, max_mutations,
                       replicates, rng, chunk_size=20000):
    """Simulating the lookup table P(Fitch score s, gene passes filter | M mutations)
    input_leaf_mask: bool list over leaves (ordered by node index), True for tips in the input
    keep_fn: function(number of input genomes with the gene) -> bool array (filter of the input
             data), or None if simulated genes are not filtered
    max_mutations: largest number of mutations M in the table
    return: m_values (array 1..max_mutations), tables {0: array, 1: array} with
            tables[r][i, s] = P(s, pass | m_values[i], root state r)
    """
    tree_arr = tree_arrays(new_tree)
    input_leaf_mask = np.asarray(input_leaf_mask, dtype=bool)
    num_input_tips = int(input_leaf_mask.sum())
    m_values = np.arange(1, max_mutations + 1)
    max_score = len(new_tree)  # Fitch score cannot exceed the number of nodes
    counts = {0: np.zeros((len(m_values), max_score + 1)),
              1: np.zeros((len(m_values), max_score + 1))}

    # process several values of M at once to keep numpy arrays reasonably large
    m_per_chunk = max(1, chunk_size // replicates)
    for start in tqdm(range(0, len(m_values), m_per_chunk)):
        chunk_m = m_values[start:start + m_per_chunk]
        mutation_sets = []
        for m in chunk_m:
            mutation_sets.extend(simulate_mutation_sets(rng, branch_dist, m, replicates))
        leaf_states, scores = tips_and_fitch(tree_arr, mutation_sets)
        num_present = leaf_states[input_leaf_mask].sum(axis=0)
        for root in (0, 1):
            # root state 1 gives the complementary tip pattern with the same Fitch score
            present = num_present if root == 0 else num_input_tips - num_present
            keep = np.ones(len(scores), dtype=bool) if keep_fn is None else keep_fn(present)
            for j, m in enumerate(chunk_m):
                block = slice(j * replicates, (j + 1) * replicates)
                counts[root][m - 1] += np.bincount(scores[block][keep[block]],
                                                   minlength=max_score + 1)

    used = max(int(np.max(np.nonzero(counts[r].sum(axis=0))[0], initial=0)) for r in (0, 1))
    tables = {r: counts[r][:, :used + 1] / replicates for r in (0, 1)}
    return m_values, tables


def estimate_mutation_distribution(fitch_counts, class_tables, max_iterations=1000,
                                   tolerance=1e-10):
    """EM estimate of the prior over M given the input Fitch score counts
    fitch_counts: {class: np.array of number of input genes per Fitch score s (index s)}
    class_tables: {class: np.array P(s, pass | M) (rows: M, columns: s)}
    return: prior (np.array over M), number of iterations, unexplained (np.array over s of input
            genes whose Fitch score never occurred in the lookup table)
    """
    num_m, num_s = next(iter(class_tables.values())).shape
    pass_prob = {c: t.sum(axis=1) for c, t in class_tables.items()}
    prior = np.full(num_m, 1.0 / num_m)
    unexplained = np.zeros(num_s)

    def step(prior):
        new = np.zeros(num_m)
        for c, table in class_tables.items():
            q = fitch_counts[c]
            joint = prior[:, None] * table
            marginal = joint.sum(axis=0)
            observed = marginal > 0
            posterior_weighted = joint[:, observed] / marginal[observed] @ q[observed]
            with np.errstate(divide="ignore", invalid="ignore"):
                new += np.where(pass_prob[c] > 0, posterior_weighted / pass_prob[c], 0)
        return new / new.sum()

    for c, table in class_tables.items():
        support = table.sum(axis=0) > 0
        unexplained += np.where(support, 0, fitch_counts[c])

    iterations = 0
    for iterations in range(1, max_iterations + 1):
        new_prior = step(prior)
        change = np.abs(new_prior - prior).max()
        prior = new_prior
        if change < tolerance:
            break
    return prior, iterations, unexplained


def calibrate(dist_one, dist_zero, dist_both, new_tree, branch_dist, input_leaf_mask, keep_fn,
              replicates, iterations, rng):
    """Distribution of the number of mutations for the simulated genes of each root class
    dist_one, dist_zero, dist_both: input Fitch score distributions {s: number of genes} for root
                                    state 1, 0 and ambiguous
    iterations: maximal number of EM iterations (1 = uniform prior)
    return: m_dist_one, m_dist_zero, m_dist_both: {M: weight}, to be used by simulation.simu and
            simulation.additional_simu_new in place of the Fitch score distributions. Weights are
            (fractional) numbers of genes to simulate before filtering.
    return: info: dict with details for the log file
    """
    input_dists = {1: dist_one, 0: dist_zero, 2: dist_both}
    max_input_score = max(max(d.keys(), default=0) for d in input_dists.values())
    num_branches = int(np.count_nonzero(branch_dist > 0))
    max_mutations = min(num_branches, max(3 * max_input_score, max_input_score + 10))

    m_values, tables = fitch_lookup_table(new_tree, branch_dist, input_leaf_mask, keep_fn,
                                          max_mutations, replicates, rng)
    # root state of the simulated genes: as reconstructed for root classes 0 and 1, and random
    # (independently for each gene) for the ambiguous class
    class_tables = {1: tables[1], 0: tables[0], 2: 0.5 * (tables[0] + tables[1])}
    num_s = tables[0].shape[1]
    fitch_counts = {}
    for c, d in input_dists.items():
        q = np.zeros(max(num_s, max_input_score + 1))
        for s, num in d.items():
            q[s] = num
        fitch_counts[c] = q
    if max_input_score + 1 > num_s:
        class_tables = {c: np.pad(t, ((0, 0), (0, max_input_score + 1 - num_s)))
                        for c, t in class_tables.items()}

    prior, used_iterations, unexplained = estimate_mutation_distribution(
        fitch_counts, class_tables, max_iterations=iterations)

    m_dists = {}
    posterior_summary = {}
    for c, table in class_tables.items():
        q = fitch_counts[c]
        pass_prob = table.sum(axis=1)
        joint = prior[:, None] * table
        marginal = joint.sum(axis=0)
        weights = np.zeros(len(m_values))
        for s in np.nonzero(q)[0]:
            if marginal[s] > 0:
                posterior = joint[:, s] / marginal[s]
                with np.errstate(divide="ignore", invalid="ignore"):
                    # genes needed before filtering to obtain the genes passing the filter
                    weights += q[s] * np.where(pass_prob > 0, posterior / pass_prob, 0)
                if s not in posterior_summary:
                    posterior_summary[s] = posterior
            else:
                # Fitch score not observed in the lookup table: fall back to M = s
                weights[min(s, len(m_values)) - 1] += q[s]
        m_dists[c] = {int(m): float(w) for m, w in zip(m_values, weights) if w > 0}

    # summary of P(M | s) for the log file (posterior of root class 0/1; root has no influence
    # without frequency filter)
    summary = {}
    for s in sorted(posterior_summary):
        posterior = posterior_summary[s]
        cdf = np.cumsum(posterior)
        summary[int(s)] = (float(posterior @ m_values),
                           int(m_values[np.searchsorted(cdf, 0.05)]),
                           int(m_values[np.searchsorted(cdf, 0.5)]),
                           int(m_values[min(np.searchsorted(cdf, 0.95), len(m_values) - 1)]))
    info = {"max_mutations": int(max_mutations), "num_branches": num_branches,
            "replicates": replicates, "iterations": used_iterations,
            "max_iterations": iterations, "unexplained_genes": int(unexplained.sum()),
            "prior_mass_at_max": float(prior[-max(1, len(prior) // 20):].sum()),
            "posterior_summary": summary}
    return m_dists[1], m_dists[0], m_dists[2], info
