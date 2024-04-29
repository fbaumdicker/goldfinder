
from matplotlib import testing

import argparse

parser = argparse.ArgumentParser(prog="goldfinder"
,description="just testing", epilog="what is this")
args = parser.parse_args()










import numpy as np

s = 60000 * 60000
#print(s)

#test = np.empty(s, dtype=np.uint64)

#test[0:2000000000] = 5

#print(test[0])



from statsmodels.distributions.empirical_distribution import ECDF

#scores = [[0,2,4,50], [2,5,6,7]]


#a = np.memmap('test.mymemmap', dtype='float32', mode='w+', shape=(1,s))
#a[0:2000000000] = 5


#cdf = ECDF(a)
#print(cdf([2,3,4,5,]))


ecdf_1 = ECDF([3, 3, 1, 4])

ecdf_2 = ECDF([3, 1, 4, 1, 1, 1])


print(ecdf_1([3, 55, 0.5, 1.5]))
print(ecdf_2([3, 55, 0.5, 1.5]))

"""
import numpy as np
import pandas as pd
import numba
from numba import njit, prange
from numba.typed import Dict

from numba import types
from numba_progress import ProgressBar
import random
from sklearn.metrics import zero_one_loss
from tqdm import tqdm

import itertools

genes = 70000
genomes = 500
nodes = 2 * genomes -1


ancestors = np.random.randint(0,2,(genes,nodes))
descendants = np.random.randint(0,2,(genes,nodes))
counts = np.random.randint(1,50,genes)


@njit(parallel=True)#, fastmath=True)
def simultaneous_score_null_dist_numba(anc, desc,counts, progress_proxy):

    print("Calculating simultaneous score for the null distribution")

    counts = np.asarray(counts)
    
    #counts = numba.typed.List(counts)
    
    ancestors = np.asarray(anc)
    descendants = np.asarray(desc)

    simu_scores = Dict.empty(key_type = types.int64, value_type= types.int64)
    for x in prange(len(ancestors)):
        simu_score = np.sum(np.multiply((ancestors[x:]-descendants[x:]),(ancestors[x]-descendants[x])),axis=1)#simultaneous_score(ancestors[x:], descendants[x:], ancestors[x],descendants[x]) 
        new_counts = counts[x:] * counts[x]
        new_counts[0] = counts[x] * (counts[x] -1) / 2
        for y in range(len(simu_score)):
            simu_scores[simu_score[y]] = simu_scores.get(simu_score[y],0) + new_counts[y]

        progress_proxy.update(1)
    return simu_scores

rows, columns = ancestors.shape

with ProgressBar(total=rows ) as progress:
        result = simultaneous_score_null_dist_numba(ancestors, descendants, counts, progress)

print(result)

@njit(parallel=True)#, fastmath=True)
def simultaneous_score_null_dist_numba_no_change(anc, desc,progress_proxy):

    print("Calculating simultaneous score for the null distribution")

    
    ancestors = np.asarray(anc)
    descendants = np.asarray(desc)

    simu_scores = Dict.empty(key_type = types.int64, value_type= types.int64)
    for x in prange(len(ancestors)):
        simu_score = np.sum(np.multiply((ancestors[x:]-descendants[x:]),(ancestors[x]-descendants[x])),axis=1)#simultaneous_score(ancestors[x:], descendants[x:], ancestors[x],descendants[x]) 
        #new_counts = counts[x:] * counts[x]
        #new_counts[0] = counts[x] * (counts[x] -1) / 2
        for y in range(len(simu_score)):
            simu_scores[simu_score[y]] = simu_scores.get(simu_score[y],0) +1 #+ new_counts[y]

        progress_proxy.update(1)
    return simu_scores

#rows, columns = ancestors.shape

#with ProgressBar(total=rows ) as progress:
#        result = simultaneous_score_null_dist_numba_no_change(ancestors, descendants,progress)

#print(result)









































def simultaneous_score(panc,pdesc,ganc,gdes):
    
    score = np.sum(np.multiply((panc-pdesc),(ganc-gdes)),axis=1)
    return score


def simultaneous_score_null_dist(anc, desc):
   

    print("Calculating simultaneous score for the null distribution")

    #counts = np.asarray(counts)
    ancestors = np.asarray(anc)
    descendants = np.asarray(desc)
    #simu_scores = {}

    simu_scores = {}
    for x in tqdm(range(len(ancestors))):
        simu_score = np.sum(np.multiply((ancestors[x:]-descendants[x:]),(ancestors[x]-descendants[x])),axis=1)#simultaneous_score(ancestors[x:], descendants[x:], ancestors[x],descendants[x]) 
        #new_counts = counts[x:] * counts[x]
        #new_counts[0] = counts[x] * (counts[x] -1) / 2
        for y in range(len(simu_score)):
            simu_scores[simu_score[y]] = simu_scores.get(simu_score[y],0) +1 # + new_counts[y]

    return simu_scores

#result = simultaneous_score_null_dist(ancestors, descendants)

"""