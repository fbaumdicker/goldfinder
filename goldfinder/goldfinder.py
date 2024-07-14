import checks
import clustering
import simulation
import output
import tree_reconstruction
import preprocessing
import data_import
import testing
import anc_recon
import scoring
import argparse
import random
import numpy as np
import tests


def argparser():
    """Store user input and return them as Namespace object"""
    parser = argparse.ArgumentParser(
        prog="goldfinder", description="Finding co-occurring genes while accounting for linkage by descent.",
        epilog = "If you use Goldfinder, please cite:"
        " Goldfinder: Unraveling Networks of Gene Co-occurrence and Avoidance in Bacterial Pangenomes;"
        #" Athina Gavriilidou, Emilian Paulitz, Christian Resl, Nadine Ziemert, Anne Kupczok, Franz Baumdicker;"
        " bioRxiv; 2024.04.29.591652; doi: 10.1101/2024.04.29.591652",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Input
    parser.add_argument("-i", "--input", help="File containing gene presence/absence data, or panX "
                        "folder", default=argparse.SUPPRESS)
    parser.add_argument("-f", "--file_type", nargs='?', choices=["roary", "tab", "panx", "matrix"],
                        default="matrix", help="Clarifying which input file type Goldfinder is "
                        "dealing with: a) roary (or panaroo), b) tab, c) panx d) .csv matrix in "
                        "pandas format")
    parser.add_argument("-t", "--tree", nargs='?', const=None,
                        help="User's phylogenetic tree in newick string format. Required, but if "
                        "given without argument, a phylogenetic tree is calculated from the input "
                        "table and used in further analysis.", default=argparse.SUPPRESS)
    # TODO test this parameter with panaroo and panx input
    parser.add_argument("-m", "--metadata", help="Path to a tab-separated table "
                        "containing metadata about the genes (e.g. functional annotation, "
                        "groupings, etc.). First column must identify the gene (i.e. coincide with "
                        "first col of input), first row must identify the type of metadata. The "
                        "metadata will be included in the results table "
                        "and occur in the visualization. If the format is roary or panaroo, columns"
                        " 'Non-unique Gene name' and 'Annotation' are used anyway.", default=None)

    # Output
    parser.add_argument("-o", "--output", nargs='?', default="./output",
                        help="Path to empty or non-existent folder where output files should be stored")
    parser.add_argument("-O", "--force_output", action="store_true", help="Set this in order to "
                        "silence the warning about the output folder existing and instead continue "
                        "execution. Warning: Might lead to problems and overwriting if files with "
                        "the same name already exist.")

    # Preprocessing
    parser.add_argument("-prep", "--preprocess", action="store_true",
                        help="Genes present in less than 5%% of sample are removed")
    parser.add_argument("-tinf", "--tree_inference", nargs='?', choices=["nj", "ml"], default="nj",
                        help="Method with which a phylogenic tree will be inferred. Choices: nj "
                        "(Neighbor joining), ml (Maximum likelihood). Default: nj")

    # Simulation
    parser.add_argument("-g", "--genes_simulated", nargs='?', default=10000,
                        help="integer specifying number of genes to be simulated for estimating the"
                        " null distribution.")
    parser.add_argument("-add", "--additional", default=0.9999, type=float, help="Percentile of "
                        "scores for which addional fine-grained simulation will be performed. Only "
                        "relevant for simultaneous score.")

    # Testing
    parser.add_argument("-s", "--score", nargs='?', choices=['terminal', 'simultaneous',
                                                             'subsequent', 'coinfinder'],
                        default='simultaneous', help="Type of score used for evaluating gene-gene "
                        "association/dissociation")
    parser.add_argument("-a", "--alpha", nargs='?', default=0.05, help="significance level")
    parser.add_argument("-pcor", "--pvalue_correction", nargs='?',
                        choices=['bonferroni', 'fdr', 'none'], default="fdr", help="p-value "
                        "correction method: bonferroni, false discovery rate, or none")
    parser.add_argument("-k", "--known_associations", nargs='?', default=0, help="Number of surely"
                        " occurring associations (or dissociations) or path to a list in the format"
                        " 'Gene1\\tGene2' (\\t denotes tab) of all surely associated gene pairs.' "
                        "Number is inferred from the list and used to improve multiple testing "
                        "correction (Bonferroni and FDR). If a list, the gene pairs will occur "
                        "in the output, even if the corresponding p-value is not siginificant.")

    # Clustering
    parser.add_argument("-inf", "--inflation", nargs='?', default=2.5,
                        help="Inflation affects granularity/resolution of clustering outcome and "
                        "needs to be >1")
    parser.add_argument("-n", "--no_clustering", action='store_true', help="Set this if no "
                        "clustering should be performed.")

    # Miscellaneous
    parser.add_argument("-c", "--coocurrence", nargs='?', choices=["association", "dissociation",
                                                                   "both"], default="association",
                        help="Whether to calculate associating or dissociating genes, or both. If "
                        "both and clustering is enabled, mean dissociation between association "
                        "clusters will be calculated and Cytoscape output will be generated.")
    parser.add_argument("--seed", nargs='?', type=int, help="An integer to be used as a random seed"
                        " in order to make the result reproducible.", default=None)
    parser.add_argument("--tests", action='store_true', help="Perform tests and exit.")

    args = parser.parse_args()

    # Check for required arguments if user does not just want to check tests
    if not args.tests:
        missing_args = []

        try:
            args.input
        except AttributeError:
            missing_args.append('-i/--input')

        try:
            args.tree
        except AttributeError:
            missing_args.append('-t/--tree')

        if missing_args:
            print('Error: the following arguments are required: ' + ', '.join(missing_args))
            exit('To display the help message, please use goldfinder.py --help')



    return args


def main():
    p = argparser()

    random.seed(p.seed)
    np.random.seed(p.seed)

    if p.tests:
        tests.perform_tests()
        exit()

    print("Checking files and parser arguments")
    checks.check_input_file(p.input, p.file_type)

    print("Creating output folder")
    output.create_output_folder(p.output, p.force_output)

    print("Loading data")
    df, locus_dict, metadata, num_known_assoc, known_assoc = data_import.load_input(
        p.input, p.file_type, p.metadata, p.known_associations)

    # Preprocessing
    pre_df = preprocessing.preproc(df, p.preprocess)

    # Constructing tree structure from nwk string or from distance matrix
    tree_struc, tip_lvlorder, branch_distances = tree_reconstruction.main(p.tree, p.tree_inference,
                                                                          pre_df, p.output)

    print("Ancestral reconstruction and determining homoplasy distribution")
    # fitch score dicts are used to determine how often which number of mutations is simulated
    (input_leaves_state, input_desc_state, input_anc_state, fitch_score_root_one,
     fitch_score_root_zero, fitch_score_root_both) = anc_recon.recon(pre_df, tree_struc,
                                                                     tip_lvlorder)

    print("Simulate for three possible root states")
    (simul_anc_states, simul_desc_states, simul_leaves, simul_counts, new_dist_one, new_dist_zero,
     new_dist_both) = simulation.simu(fitch_score_root_one, fitch_score_root_zero,
                                      fitch_score_root_both, tree_struc, branch_distances,
                                      int(p.genes_simulated))

    # Calculate scores for simulated data and input data")
    null_dist_scores, input_scores = scoring.scoring_procedure(
        p.score, simul_anc_states, simul_desc_states, simul_leaves, simul_counts,
        fitch_score_root_one, fitch_score_root_zero, fitch_score_root_both, tree_struc,
        branch_distances, int(p.genes_simulated), pre_df, input_anc_state, input_desc_state,
        input_leaves_state, p.output)

    # Additional simulation of distribution's tail so far only implemented for simulatenous score
    if p.score == "simultaneous":
        print("Additional Simulation for extreme region(s)")
        (add_simul_anc_states, add_simul_desc_states, add_simul_leaves, add_simul_counts, k,
         k_frac) = simulation.additional_simu_new(
             null_dist_scores, fitch_score_root_one, fitch_score_root_zero, fitch_score_root_both,
             tree_struc, branch_distances, int(p.genes_simulated), p.output, p.additional)

        if k_frac == 1.0:
            print("Warning: No additional simulation performed. Please use a larger percentile.")
        else:
            print("Scoring of additional simulation")
            null_dist_scores = scoring.additional_scoring(
                add_simul_anc_states, add_simul_desc_states, add_simul_counts, k, k_frac,
                null_dist_scores, p.output)

    cluster_dict = None
    modes = ["association", "dissociation"] if p.coocurrence == 'both' else [p.coocurrence]
    for mode in modes:
        print(f"Testing for {mode}")
        significant_score_indices, p_values_unadj, p_values_adj, sig_lvl, known_assoc = \
            testing.testing_procedure(null_dist_scores, input_scores, mode, p.alpha,
                                      p.pvalue_correction, p.output, num_known_assoc, known_assoc)

        if p.coocurrence == 'both' and mode == 'dissociation' and not p.no_clustering:
            print("Calculating Average Dissociation between Clusters")
            dissoc_freq = clustering.dissociation_freq(cluster_dict, p_values_adj,
                                                       p.file_type in ["tab", "matrix"])

            # Write to file, along with global fraction of significantly dissociated gene pairs
            num_significant = significant_score_indices[0].size
            num_gene_pairs = input_scores.shape[0] * (input_scores.shape[0] - 1) / 2
            output.cluster_dissoc(dissoc_freq, num_significant / num_gene_pairs, p.output)

            output.clusters_cytoscape(dissoc_freq, cluster_dict, p.output)

        perform_clustering = mode == "association" and not p.no_clustering
        if perform_clustering:
            print("Clustering")
        cluster_dict, clusters = clustering.cluster_procedure(p_values_adj, sig_lvl,
                                                              float(p.inflation), p.output,
                                                              perform_clustering)

        print(f"Preparing output files for {mode}")
        output.result_procedure(p_values_adj, p_values_unadj, significant_score_indices,
                                cluster_dict, clusters, locus_dict, p.output, p.score, mode,
                                p.file_type, perform_clustering, metadata,
                                known_assoc, p.coocurrence == 'both' and not p.no_clustering)

    print("Analysis is finished")


if __name__ == "__main__":
    main()
