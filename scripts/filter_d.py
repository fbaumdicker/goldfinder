import pandas as pd

wd = './coin_time/'
threshold = -0.2

D = pd.read_csv(f'{wd}coin_nodes.tsv', sep='\t').set_index('ID')
Ddict = D['Result'].to_dict()

df = pd.read_csv(f'{wd}coin_pairs.tsv', sep='\t')

# Filter out pairs where either gene has no D-value or either gene's D-value is below the threshold
df = df.loc[df['Source'].isin(D.index) & df['Target'].isin(D.index), :].copy()

# Calculate D value for gene pairs
df['sourceD'] = df['Source'].map(lambda x: Ddict[x])
df['targetD'] = df['Target'].map(lambda x: Ddict[x])
df['minD'] = df[['sourceD', 'targetD']].min(axis=1)

df.loc[df['minD'] > threshold, :].drop(['sourceD', 'targetD', 'minD'], axis=1).to_csv(f'{wd}coin_filtered_d.tsv', sep='\t', index=False)
