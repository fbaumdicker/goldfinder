import markov_clustering as mc
import pandas as pd
import numpy as np

wd = 'C:/Users/emilp/Documents/Uni/Forschungsprojekt/'
in_file = f'{wd}data/time_results/coin_filtered_d_p.tsv'
out_file = f'{wd}data/time_results/coin_clusters_new.txt'

df = pd.read_csv(in_file, sep='\t')[['Source', 'Target', 'p']]

# get all genes, and a dict gene: idx
genes = list(df[['Source', 'Target']].stack().unique())
dim = len(genes)
gene_toi = {gene: i for i, gene in enumerate(genes)}

# translate gene names to their idx
df['Source_i'] = df['Source'].apply(lambda gene: gene_toi[gene])
df['Target_i'] = df['Target'].apply(lambda gene: gene_toi[gene])

print('Create and fill adjacency matrix')
adj_mtx = np.array([[0.0] * dim] * dim)
for i, row in df.iterrows():
    adj_mtx[row['Source_i']][row['Target_i']] = row['p']
    adj_mtx[row['Target_i']][row['Source_i']] = row['p']

# set distance of all significant pairs to 1, all other pairs to 0
(adj_mtx != 0).astype(int)

print('Run MCL')
result = mc.run_mcl(adj_mtx, inflation=1.3)
clusters = mc.get_clusters(result)

# Filter results: No duplicate clusters
st = set()
filtered_clusters = []
for tup in clusters:
    filtered_tup = tuple(gene for gene in tup if gene not in st)
    st.update(filtered_tup)
    if len(filtered_tup) > 1:
        filtered_clusters.append(filtered_tup)

clusters = filtered_clusters

print('Output in same format as Goldfinder')
with open(out_file, 'w') as f:
    gene_cluster_nr = 0
    for x in range(len(clusters)):
        if len(clusters[x]) > 1:
            cluster_size = len(clusters[x])
            gene_cluster_nr += 1
            s = ">" + str(gene_cluster_nr) + "," + str(cluster_size) + "\n"
            f.write(s)
            for gene_loc_i in range(cluster_size):
                s = genes[clusters[x][gene_loc_i]]

                if gene_loc_i < cluster_size - 1:  # not the last gene
                    s += ','

                f.write(s + '\n')
