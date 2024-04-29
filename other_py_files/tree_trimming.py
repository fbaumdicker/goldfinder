
import matplotlib.pyplot as plt
import pandas as pd
from ete3 import Tree

"""
cluster_size = []

with open("/home/chris/Workspace/coin_treewas_mp/output/clusters.txt") as f:

    lines = f.readlines()


    for line in lines:
        if line[0] == ">":
            cluster_size.append(line.split(",")[1][0])
    f.close()

cluster_size.sort()
print(cluster_size)
plt.hist(cluster_size, bins = 10)
plt.xlabel("Cluster size")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig("cluster_size_hist.png")

"""

df = pd.read_csv("/home/chris/Downloads/evidence_of_selection_pseudomonas/gene_presence_absence.csv", index_col=0, low_memory=False)

annotation = df.iloc[:,0:2]
df = df.drop(df.iloc[:, 0:13], axis = 1) #remaining columns are samples/species
df = df.notnull().astype('int') #empty cells'value replaced with 0 and non-empty with 1
    
#incorporating annotation in gene_name
tmp = {}
for index, row in annotation.iterrows():
    tmp[index] = str(index)+"/"+str(row[0])+"/"+str(row[1]).replace(",", "")
#renaming gene name to gene name + annotation
df.rename(index = tmp, inplace=True)

genomes = list(df.columns)
genomes_dict = dict.fromkeys(genomes, 0)

nwk = open("/home/chris/Downloads/evidence_of_selection_pseudomonas/core.gene.align.fasta-gb.treefile").read().replace('t', 'V')
t = Tree(nwk)
tdict = {}
for leaf in t.iter_leaf_names():
    tdict[leaf] = 0


tdict_not_in_genomes_dict = 0

for key in tdict.keys():
    if key not in genomes_dict:
        tdict_not_in_genomes_dict += 1

print("tdict size: ", len(tdict), " Not in genomes_dict: ", tdict_not_in_genomes_dict)

genomes_dict_not_in_tdict = 0

for key in genomes_dict.keys():
    if key not in tdict:
        genomes_dict_not_in_tdict += 1

print("genomes_dict size: ", len(genomes_dict), " Not in tdict: ", genomes_dict_not_in_tdict)