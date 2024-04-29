import argparse
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

import scoring
import anc_recon
import testing
import input
import preprocessing
import tree_reconstruction
import output
import simulation
import assoc_net
import checks


def argparser ():
    """Store user input and return them as Namespace object"""
    parser = argparse.ArgumentParser(prog = "goldfinder", description="Finding co-occurring genes while accounting for linkage by descent.")
    parser.add_argument("-i", "--input", required=True, help="File containing gene presence/absence data")
    parser.add_argument("-f", "--file_type", nargs='?', choices=["roary", "tab", "panx", "panaroo"], default="roary", help="Clarifying which input file type Goldfinder has to deal with")
    parser.add_argument("-p", "--preprocess", action="store_true", help="Genes present in less than 5%% of sample are removed")
    parser.add_argument("-t", "--tree",required=True, nargs= '?', const=None, help="User's phylogenetic tree in newick string format")
    parser.add_argument("-g", "--genes_simulated", nargs='?', default = 10000, help ="integer specifying number of genes to be simulated for estimating the null distribution" )
    parser.add_argument("-a", "--alpha", nargs='?', default=0.05, help = "significance level")
    parser.add_argument("-pcor", "--pvalue_correction", nargs='?', choices = ['bonferroni', 'fdr', ""],default="fdr", help="p-value correction method: bonferroni or false discovery rate")
    parser.add_argument("-c", "--coocurrence", nargs='?', choices =['association', "dissociation"], default="association", help="String indicating either association or dissociation")
    parser.add_argument("-tinf","--tree_inference", nargs='?',choices =['nj', "ml"], default= 'nj', help="Method(Neighbor joining, Maximum likelihood) with which a phylogenic will be inferred")
    parser.add_argument("-o", "--output", nargs='?', default="./output", help="path to where output should be stored")
    parser.add_argument("-s", "--score", nargs='?', choices = ['terminal', 'simultaneous', 'subsequent', 'coinfinder'], default= 'simultaneous', help = "Type of score used for evaluating gene-gene association/dissociation" )
    parser.add_argument("-inf", "--inflation", nargs='?', default = 2, help = "Inflation affects granularity/resolution of clustering outcome")
    parser.add_argument("-n", "--network", nargs='?', choices=["yes", "no"], default="yes", help = "If an associaton network should be built in case association is measured" )

    args = parser.parse_args()
    return args

def main():
    p = argparser()

    print("Checking files and parser arguments")
    checks.check_input_file(p.input, p.file_type, p.tree)

    print("Creating output folder")
    output.create_output_folder(p.output)

    print("Loading data")
    df, locus_dict = input.load_input(p.input, p.file_type)

    print("Preprocessing")
    pre_df = preprocessing.preproc(df, p.preprocess)

    print("Constructing tree structure")
    tree_struc, tip_lvlorder, branch_distances = tree_reconstruction.tree_recon(p.tree, p.tree_inference, pre_df, p.output)

    print("Ancestral reconstruction and determining homoplasy distribution")
    input_leaves_state, input_desc_state, input_anc_state, fitch_score_root_one, fitch_score_root_zero, fitch_score_root_both = anc_recon.recon(pre_df, tree_struc, tip_lvlorder)

    print("Simulation")
    simul_anc_states, simul_desc_states, simul_leaves, simul_counts, new_dist_one, new_dist_zero, new_dist_both = simulation.simu(fitch_score_root_one, fitch_score_root_zero, fitch_score_root_both, tree_struc, branch_distances, int(p.genes_simulated))    
    
    print("self simulated dataframe")
    new_simul_leaves = []
    for x in range(len(simul_counts)):
        for y in range(simul_counts[x]):
            new_simul_leaves.append(simul_leaves[x])


    pre_dataframe = pd.DataFrame(new_simul_leaves, columns = list(pre_df.columns), index=["g" + str(x) + "/whatever/whatever" for x in range(len(new_simul_leaves))])    

    print("Preprocessing")
    pre_df = preprocessing.preproc(pre_dataframe, p.preprocess)

    print("Constructing tree structure")
    tree_struc, tip_lvlorder, branch_distances = tree_reconstruction.tree_recon(p.tree, p.tree_inference, pre_df, p.output)

    print("Ancestral reconstruction and determining homoplasy distribution")
    input_leaves_state, input_desc_state, input_anc_state, fitch_score_root_one, fitch_score_root_zero, fitch_score_root_both = anc_recon.recon(pre_df, tree_struc, tip_lvlorder)

    print("Simulation")
    simul_anc_states, simul_desc_states, simul_leaves, simul_counts, new_dist_one, new_dist_zero, new_dist_both = simulation.simu(fitch_score_root_one, fitch_score_root_zero, fitch_score_root_both, tree_struc, branch_distances, int(p.genes_simulated))    
    

    print("Scoring")
    null_dist_scores, input_scores = scoring.scoring_procedure(p.score, simul_anc_states, simul_desc_states, simul_leaves, simul_counts, fitch_score_root_one, fitch_score_root_zero, fitch_score_root_both, tree_struc, branch_distances, int(p.genes_simulated), pre_df, input_anc_state, input_desc_state, input_leaves_state, p.output)

    print("Testing")
    significant_score_indices, p_values_unadj, p_values_adj, sig_lvl = testing.testing_procedure(null_dist_scores, input_scores, p.coocurrence, p.alpha, p.pvalue_correction, p.output)


    print("Clustering")
    cluster_dict, clusters = assoc_net.cluster_procedure(p_values_adj, sig_lvl, float(p.inflation), p.output, p.file_type, p.network)

    print("Preparing output files")
    output.result_procedure (p_values_adj, p_values_unadj, significant_score_indices, cluster_dict, clusters, locus_dict, p.output, p.score, p.file_type, p.network)

    print("Analysis is finished")


if __name__ == "__main__":
    main()
