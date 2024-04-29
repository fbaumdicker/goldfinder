import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.stats import pearsonr
from tqdm import tqdm

import networkx as nx
import markov_clustering as mc
import matplotlib.pyplot as plt

#coinfinder main paper input
#cf_d_all = pd.read_csv("~/Workspace/coincident_nodes.tsv",index_col=0, sep="\t")
#evidence for selection paper input
cf_d_all = pd.read_csv("~/Downloads/evidence_of_selection_pseudomonas/pseudomonas-manuscript-main/all-abundant-accessory-genes_nodes.tsv",index_col=0, sep="\t")
cf_d_all_2 = pd.read_csv("~/Downloads/evidence_of_selection_pseudomonas/pseudomonas-manuscript-main/coinfinder_L0.05cmgfspGTRIGiU0.9F0.05_D04_nodes.tsv",index_col=0, sep="\t")
#making d value dictionary with all data set
cf_d_first = cf_d_all["Result"].tolist()
cf_d_second = cf_d_all.index.tolist()
cf_d_first_2 = cf_d_all_2["Result"].tolist()
cf_d_second_2 = cf_d_all_2.index.tolist()
cf_d = {}
for x in range(len(cf_d_first)):
    cf_d[cf_d_second[x]] = cf_d_first[x]
for x in range(len(cf_d_first_2)):
    cf_d[cf_d_second_2[x]] = cf_d_first_2[x]


#goldfinder output for coinfinder main paper input
#gf_fdr = pd.read_csv("./output/simultaneous_score_significant_gene_pairs.txt", header = None, sep = ",")
#goldfinder output for evidecne for selection paper input
gf_fdr = pd.read_csv("~/Downloads/evidence_of_selection_pseudomonas/goldfinder_output_evidence_for_selection/simultaneous_significant_gene_pairs.txt", header = None, sep = ",")
#gf_fdr = pd.read_csv("./output/simultaneous_significant_gene_pairs.txt", header = None, sep = ",")
gf_first_fdr = list(gf_fdr[0])
#gf_second_fdr = list(gf_fdr[1])
gf_second_fdr = list(gf_fdr[3])
#gf_p_values_fdr = list(gf_fdr[2])
gf_p_values_fdr = list(gf_fdr[7])

#unique gene names for using as column index names
gf_joined_fdr = gf_first_fdr + gf_second_fdr
gf_joined_fdr = list(set(gf_joined_fdr))

#print(gf_joined_fdr)

count = 0
for x in range(len(gf_first_fdr)):
    if gf_first_fdr[x] in cf_d and gf_second_fdr[x] in cf_d:
        if cf_d[gf_first_fdr[x]] >= -0.4 and cf_d[gf_second_fdr[x]] >= -0.4:
            count += 1
print("Count: ", count)


gf_d_values_fdr = []
gf_d_values_greater_05 = []
for x in range(len(gf_joined_fdr)):
    gf_d_values_fdr.append(cf_d.get(gf_joined_fdr[x],100))
    if cf_d.get(gf_joined_fdr[x],0)>= -0.4:
        gf_d_values_greater_05.append(gf_joined_fdr[x])
    
#print("how many fdr genes of Goldfinder not in d value list: ", gf_d_values_fdr.count(100))
gf_d_values_fdr = [x for x in gf_d_values_fdr if x != 100]
plt.hist(gf_d_values_fdr, bins= 20, alpha= 0.2, label="goldfinder")
plt.savefig("d_values_goldfinder_fdr")
plt.close()


#dictionary to store efficiently p values and names

#adjacency matrix
#ad = pd.DataFrame(10, index = gf_joined_fdr, columns=gf_joined_fdr)
#for x in tqdm(range(len(gf_first_fdr))):
#    ad.loc[gf_first_fdr[x], gf_second_fdr[x]] = gf_p_values_fdr[x]
#    ad.loc[gf_second_fdr[x], gf_first_fdr[x]] = gf_p_values_fdr[x]

#ad[ad != 10] = 1
#ad[ad == 10] = 0

#ad_np = ad.to_numpy()
"""
result = mc.run_mcl(ad_np, inflation=1.4)
clusters = mc.get_clusters(result) 



with open("newtork_clusters.txt", 'a') as f:
    for cluster in tqdm(clusters):
        s=[]
        for ele in cluster: 
            s.append(gf_joined_fdr[ele])        
        e = "-".join(map(str,s)) + "\n"
        f.write(e)
    f.close()


mc.draw_graph(ad_np, clusters, with_labels=False)
plt.savefig("network_goldfinder_assoc")
plt.close()

cluster_size = {}
for x in tqdm(clusters):
    cluster_size[len(x)] = cluster_size.get(len(x),0)+1
print(cluster_size)

#plt.bar(list(cluster_size.keys()), cluster_size.values())

cd_list = [key for key, val in cluster_size.items() for _ in range(val)]
plt.hist(cd_list, bins = 20)
plt.xlabel('Cluster size', fontsize=10)
plt.ylabel('Number of clusters', fontsize=10)
plt.savefig("network_goldfinder")
plt.close()
"""


#print(ad.loc["soxS", "ykoD_3"])
#print(ad["soxS", "ykoD_3"])

#soxS,ykoD_3,0.04809162467670225


#######################################Coinfinder clusters
#original coinfinder data input
#cf_pairs_bonf = pd.read_csv("~/Workspace/d_value_test/recent_coinfinder_run/coincident_pairs.tsv", sep="\t")
#evidence for selection paper input
cf_pairs_bonf = pd.read_csv("~/Downloads/evidence_of_selection_pseudomonas/coincident_pairs.tsv", sep="\t")

#unique gene names for using as column index names
cf_first = cf_pairs_bonf["Source"].tolist()
cf_second = cf_pairs_bonf["Target"].tolist()
cf_p_values = cf_pairs_bonf["p"].tolist()

count = 0
for x in range(len(cf_first)):
    if cf_first[x] in cf_d and cf_second[x] in cf_d and cf_p_values[x] <= 3.25 * pow(10, -6):
        if cf_d[cf_first[x]] >= -0.4 and cf_d[cf_second[x]] >= -0.4:
            count += 1
print("Count: ", count)

cf_joined = cf_first + cf_second
cf_joined = list(set(cf_joined))

cf_d_values = []
cf_d_values_greater_05 = []
for x in range(len(cf_joined)):
    cf_d_values.append(cf_d.get(cf_joined[x],100))
    #original coinfinder data
    #if cf_d.get(cf_joined[x],0) >= 0.5:
    #evidence of paper datat

    if cf_d.get(cf_joined[x],0) >= -0.4:
        cf_d_values_greater_05.append(cf_joined[x])
    

cf_d_values = [x for x in cf_d_values if x != 100]
plt.hist(cf_d_values, bins= 20, alpha= 0.2, label="coinfinder")
plt.savefig("d_values_coinfinder")
plt.close()


plt.hist(cf_d_values, bins= 30, alpha= 0.2, label="coinfinder")
plt.hist(gf_d_values_fdr, bins= 30, alpha= 0.2, label="goldfinder")
plt.xlabel('D-value', fontsize=10)
plt.ylabel('Number of genes', fontsize=10)
plt.legend(loc='upper right')
plt.savefig("d_values_frequency_gold_coin")
plt.close()


#number of genes in goldfinder above D-value of 0.5 that Coinfinder does not find

print("Number of genes in Goldfinder above D-value 0.5: ", len(gf_d_values_greater_05))
print("Number of genes in Coinfinder above D-value 0.5: ", len(cf_d_values_greater_05))

gf_not_in_cf = 0
for x in gf_d_values_greater_05:
    if x not in cf_d_values_greater_05:
        gf_not_in_cf += 1
print("Genes above 0.5 D in Goldfinder but not in Coinfinder: ", gf_not_in_cf)

#print(gf_d_values_greater_05)
#print(" ")
#print(cf_d_values_greater_05)

"""
cf_pairs = pd.read_csv("~/coincident_pairs.tsv", sep="\t")

#making d value dictionary with all data set
cf_first = cf_pairs["Source"].tolist()
cf_second = cf_pairs["Target"].tolist()
cf_p_values = cf_pairs["p"].tolist()

#unique gene names for using as column index names
cf_joined = cf_first + cf_second
cf_joined = list(set(cf_joined))


#adjacency matrix
cd = pd.DataFrame(10, index = cf_joined, columns=cf_joined)
for x in tqdm(range(len(cf_first))):
    cd.loc[cf_first[x], cf_second[x]] = cf_p_values[x]
    cd.loc[cf_second[x], cf_first[x]] = cf_p_values[x]

cd[cd != 10] = 1
cd[cd == 10] = 0

cd_np = cd.to_numpy()

cd_result = mc.run_mcl(cd_np, inflation=1.4)
cd_clusters = mc.get_clusters(cd_result) 
mc.draw_graph(cd_np, cd_clusters, with_labels=False)
plt.savefig("network_coinfinder_assoc")
plt.close()

cd_cluster_size = {}
for x in tqdm(cd_clusters):
    cd_cluster_size[len(x)] = cd_cluster_size.get(len(x),0)+1
print("")
print(cd_cluster_size)

cd_list = [key for key, val in cd_cluster_size.items() for _ in range(val)]
plt.hist(cd_list, bins = 20)
#plt.bar(list(cd_cluster_size.keys()), cd_cluster_size.values())


plt.xlabel('Cluster size', fontsize=10)
plt.ylabel('Number of clusters', fontsize=10)
plt.savefig("network_coinfinder")
plt.close()
"""


unique_to_gf = []
unique_to_cf = []
in_both = []


for x in tqdm(range(len(gf_joined_fdr))):
    if gf_joined_fdr[x] in cf_joined:
        in_both.append(cf_d.get(gf_joined_fdr[x],100))
    else:
        unique_to_gf.append(cf_d.get(gf_joined_fdr[x],100))

for y in tqdm(range(len(cf_joined))):
    if cf_joined[y] not in gf_joined_fdr:
        unique_to_cf.append(cf_d.get(cf_joined[y],100))

unique_to_cf = [x for x in unique_to_cf if x != 100]
unique_to_gf = [x for x in unique_to_gf if x != 100]
in_both = [x for x in in_both if x != 100]

print("unique to goldfinder: ", len(unique_to_gf))
print(" ")
print("found in goldfinder and coinfinder: ", len(in_both))
print(" ")
print("unqiue to coinfinder: ", len(unique_to_cf))
    


print("")   
print("number of unique genes in coinfinder: ", len(cf_joined))
print("number of unqiue genes in goldfinder: ", len(gf_joined_fdr))

plt.figure()
plt.hist([unique_to_cf,unique_to_gf, in_both],color = ["red", "blue", "violet"], bins= 20, stacked=True)
plt.legend({"coinfinder": "red", "goldfinder": "blue", "in both": "violet"})
#plt.xlim(xmin=-1, xmax=1)
plt.xlabel("D-value")
plt.ylabel("Number of genes")
plt.savefig("stacked_d_values")
plt.close()