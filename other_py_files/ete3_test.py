#from ete3 import Tree
import numpy as np
import random
from itertools import repeat
from tqdm import tqdm




#phylogenetic tree inference
from skbio import DistanceMatrix
from skbio.tree import nj

ids = list('abcde')
data = [[0,5,9,9,8],
        [5,0,10,10,9],
        [9,10,0,8,7],
        [9,10,8,0,3],
        [8,9,7,3,0]]

dm = DistanceMatrix(data, ids)
print(dm)

tree = nj(dm)
print(tree.ascii_art())

newick_str = nj(dm, result_constructor=str)
print(newick_str)

new_tree = tree.root_at_midpoint()
print(new_tree.ascii_art())
print(new_tree)
nts = str(new_tree)
print(nts.replace("root", ""))

    
















"""
#Association test
import pandas as pd

data = pd.read_csv("/home/chris/Workspace/coin_treewas_mp/terminal_score_significant_gene_pairs.txt", header = True)
print(data)
#data.drop(index=data.index[0], axis=0, inplace=True)
print(list(data.columns()))

"""


"""


nwk = "((PagglomeransE325:0.000000, (Epyrifoliae:0.000000, (Etasmaniensis:0.000000, PvagansC91:-0.000000):0.000000):0.232616):0.006648, ((PananatisLMG:0.000000, PantoeaspaBvalens:0.232616):0.013297, (EbillingiaeEb661:0.250345, PantoeaspAt9b:0.000000):0.013297):0.006648, EamylovoraCFBP1430:0.000000);"


test = Tree(nwk)

lvl = 0
for node in test.traverse():
    node.add_feature("gene_presence_absence", 1)
    node.add_feature("lvl_order", lvl)
    lvl += 1
    anc = node.get_ancestors()
    print(node.name, "its ancestors:", anc)
    print("------------------------")
    #print(node.gene_presence_absence)


for node in test.traverse("postorder"):
    children = node.get_descendants()
    #if len(children)>0:
    if not node.is_leaf():
        intersection = list(set([children[0].gene_presence_absence, children[1].gene_presence_absence]))
        node.add_feature("gene_presence_absence", intersection)
    print(node.name, node.gene_presence_absence)
"""  

"""
for node in test.traverse("postorder"):
    print(node.name)
    children = node.get_descendants()
    print("its children:", children)
    if not children:
        print(children[0])
    print("-------------")
"""

#s = 10000

#print(test.iter_leaves())
#[print(node.name, node.dist) for node in test.iter_leaves()]

    
"""
for x in tqdm(range(s)):
        root = random.choice([0,1]) #determining root state
        
        #tm = Tree(nwk)
        lvlnum=0 #variable for labeling tree according to level order traversal
        
        for node in test.traverse():    
            #if the node was not yet visited, add features
            if not hasattr(node,"level_order_num"):
                node.add_feature("level_order_num", lvlnum)
            if not hasattr(node,"gene_presence"):
                node.add_feature("gene_presence", root)
                lvlnum += 1
"""                 
    
#node2leaves = test.get_cached_content()
#print(test)
#for n in test.traverse():
#    print(n.name, node2leaves[n])

#t = Tree(nwk) #an ete3 object of the phylogenetic tree
#t2 = Tree(nwk)

#print(t2.get_leaves())

"""
root = random.choice([0,1])
t2.add_feature("gene_presence", root)
#print(t2.gene_presence)
for node in t2.traverse():
    node.add_feature("gene_presence", 10)
print(t2.gene_presence)

descendants = t2.get_descendants()

[x.add_feature("gene_presence",root) if not hasattr(x, "gene_presence") else x.add_feature("gene_presence", 20) for x in descendants]



for node in t2.traverse():
    print(node.gene_presence)

"""
"""
for node in t.traverse():
    print(node.name)
    print("node feature", node.features)
    node.add_feature("gene_presence",1)
    print("new node feature", node.features)
    print(node.gene_presence)
"""


"""
print("----------------------")
for node in t.traverse():
    print(node.gene_presence)

print("-------------------")

for node in t2.traverse():
    node.gene_presence = 2
    print(node.gene_presence)

for node in t.traverse():
    print(node.gene_presence)
"""

#for node in t.traverse():
#    print(node.name, ": ", node.children)

#trees = [Tree(nwk)]*10000

#for tree in trees:
#    for node in tree.traverse():
#        node.add_feature("gene_presence", 1)

#for node in trees[0]:
#    #print(node.gene_presence)
#    t = node.get_children()
#    print(t)

#get_descendants()
#iter_descendants()
"""
for node in t2.traverse():
    n = node.get_ancestors()
    print(n)
    for c in n:
        for child in c.traverse():
            print(child.name)
        break
    print("--------------")
"""
"""
num = 0
for node in t2.traverse():
    node.add_feature("level_order_numb", num)
    num += 1
    print(node.name, node.level_order_numb)
"""
#print(np.arange(8))