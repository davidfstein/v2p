import pandas as pd
import sys
import os

on_the_fly = pd.DataFrame()
precomputed = pd.DataFrame()

if os.path.exists(sys.argv[1]):
    precomputed = pd.read_csv(sys.argv[1])
if os.path.exists(sys.argv[2]):
    on_the_fly = pd.read_csv(sys.argv[2])

out = pd.concat([precomputed, on_the_fly], ignore_index=True)
out['CHROM'] = out['ID'].apply(lambda v: v.split('_')[0])
out['POS'] = out['ID'].apply(lambda v: v.split('_')[1])
out = out.sort_values(['CHROM', 'POS'])
out = out.drop(columns=['CHROM','POS'])
out = out[['ID'] + out.columns.tolist()[:-1]]
out.to_csv(sys.argv[2], index=None)
