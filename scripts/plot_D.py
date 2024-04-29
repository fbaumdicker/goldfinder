import pandas as pd
import matplotlib.pyplot as plt

wd = 'coin_time/'
D = pd.read_csv(f'{wd}coin_nodes.tsv', sep='\t').set_index('ID')
Ddict = D['Result'].to_dict()

df = pd.read_csv(f'{wd}coin_pairs.tsv', sep='\t')
# Filter out pairs where either gene has no D-value
df = df.loc[df['Source'].isin(D.index) & df['Target'].isin(D.index), :].copy()

# Histogram
genes = df[['Source', 'Target']].stack().unique()
(pd.Series(genes)
   .map(lambda x: Ddict[x])
   .loc[lambda x: (x > -1) & (x < 0.5)]
   .hist(bins=100, grid=False, color='dimgray', edgecolor='white'))
plt.savefig(f'{wd}D_hist.png', dpi = 200)
plt.clf()

# Calculate D value for gene pairs
df['sourceD'] = df['Source'].map(lambda x: Ddict[x])
df['targetD'] = df['Target'].map(lambda x: Ddict[x])
df['minD'] = df[['sourceD', 'targetD']].min(axis=1)

# Gene Pair Plot
x = [i/10 - 1 for i in range(21)]
y = [df.loc[df['minD'] > threshold].shape[0] for threshold in x]

plt.plot(x, y, 'k.')
#plt.vlines([-0.3], plt.ylim()[0], plt.ylim()[1], linestyles='dashed', color = 'k')
plt.xticks(x[::2])
plt.xlabel('D-value')
plt.ylabel('# of coincident pairs')
plt.savefig(f'{wd}D_scatter.png', dpi = 200)
