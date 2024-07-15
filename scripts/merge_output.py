import pandas as pd
import sys

precomputed = pd.read_csv(sys.argv[1])
on_the_fly = pd.read_csv(sys.argv[2])
on_the_fly = on_the_fly[precomputed.columns]
out = pd.concat([precomputed, on_the_fly], ignore_index=True)
out['CHROM'] = out['ID'].apply(lambda v: v.split('_')[0])
out['POS'] = out['ID'].apply(lambda v: v.split('_')[1])
out = out.sort_values(['CHROM', 'POS'])
out = out.drop(columns=['CHROM','POS'])
out.to_csv(sys.argv[2], index=None)
