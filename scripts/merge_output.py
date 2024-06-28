import pandas as pd
import sys

on_the_fly = pd.read_csv(sys.argv[2])
precomputed = pd.read_csv(sys.argv[1])

out = pd.concat([precomputed, on_the_fly])
out['CHROM'] = out['ID'].apply(lambda v: v.split('_')[0])
out['POS'] = out['ID'].apply(lambda v: v.split('_')[1])
out = out.sort_values(['CHROM', 'POS'])
out = out.drop(columns=['CHROM','POS'])
out = out[['ID'] + out.columns.tolist()[:-1]]
out.to_csv(sys.argv[2], index=None)