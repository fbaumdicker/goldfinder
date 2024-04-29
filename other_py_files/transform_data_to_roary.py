import pandas as pd
import numpy as np

df = pd.read_csv("/home/chris/Workspace/coin_treewas_mp/example_files/simulated_data/shuffled_presence_ecolitree.csv", index_col=0, low_memory=False)
df.replace(0, np.nan, inplace=True)
df.replace(1, "x", inplace=True)
print(df)
df.to_csv("shuffled_updated.csv", sep=",")