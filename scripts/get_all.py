import sys
import pandas as pd
from get_proteinID_features import get_proteinID_features
from get_geneID_features import get_geneID_features
from get_aa_features import get_aa_features
from get_cadd_features import get_cadd_features

def get_all(in_path, out_path, n_cores=1):
    proteinID_features = get_proteinID_features(in_path, n_cores)
    geneID_features = get_geneID_features()
    aa_features = get_aa_features(in_path)
    cadd_features = get_cadd_features(in_path, n_cores)
    assert(proteinID_features.shape[0] == geneID_features.shape[0] == aa_features.shape[0])
    out = pd.concat([proteinID_features, geneID_features, aa_features], axis=1)
    out['ID'] = pd.read_parquet(in_path, columns=['ID'])
    out = out.merge(cadd_features, on='ID', how='left')
    out['MMSp_donorIntron'] = out['MMSp_donorIntron'].replace('-', None)
    out = out.replace('NULL', None)
    out = out.replace('.', None)
    out = out.drop(columns=['Chrom'])
    try:
        out.to_parquet(out_path, index=None)
    except Exception as e:
        print(e)
        print('Warning: Could not save as parquet, falling back to csv.')
        out.to_csv(out_path, index=None)
    return out

if __name__  == '__main__':
    get_all(sys.argv[1], sys.argv[2], sys.argv[3])
