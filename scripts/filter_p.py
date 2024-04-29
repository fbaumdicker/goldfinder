import pandas as pd

wd = 'coin_time/'
p_thresh = 0.05

df = pd.read_csv(f'{wd}coin_filtered_d.tsv', sep='\t')
df.loc[df['p'] < 0.05, :].to_csv(f'{wd}coin_filtered_d_p.tsv', sep='\t', index=False)
