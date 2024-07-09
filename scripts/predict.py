import pandas as pd
import pyarrow.parquet as pq
import numpy as np
from tqdm import tqdm
import multiprocessing as mp
from glob import glob
import sys
import os
import itertools
import joblib

PREFIX = 'full_'
POSTFIX = '_full'
PROJECT_DIR = os.environ['V2P_DIR']
OUT_COLUMNS = ['HP:0033127','HP:0040064','HP:0000707','HP:0001939','HP:0000152','HP:0001626','HP:0000119','HP:0000478',
                'HP:0002715','HP:0001574','HP:0001871','HP:0025031','HP:0002664','HP:0002086','HP:0000818','HP:0000598',
                'HP:0025354','HP:0001197','HP:0001507','HP:0025142','HP:0000769','HP:0001608','HP:0045027','Pathogenic']

CUTOFFS = {'HP:0033127': 0.13375073118275554, 'HP:0040064': 0.1027121652336587, 'HP:0000707': 0.1280226940316797, 'HP:0001939': 0.13289699594521331, 'HP:0000152': 0.10556380274414325, 'HP:0001626': 0.1553588637596261, 'HP:0000119': 0.1332160134368378, 'HP:0000478': 0.20826082412161753, 'HP:0002715': 0.15873701583024458, 'HP:0001574': 0.16293394524364474, 'HP:0001871': 0.13980646528040241, 'HP:0025031': 0.1296962674378072, 'HP:0002664': 0.10575180529819386, 'HP:0002086': 0.02704949652042875, 'HP:0000818': 0.12326556052836206, 'HP:0000598': 0.12354089044413866, 'HP:0025354': 0.23111810436978011, 'HP:0001197': 0.0020001575650973, 'HP:0001507': 0.21537008781720823, 'HP:0025142': 0.0988696186665396, 'HP:0000769': 0.13120347198135124, 'HP:0001608': 0.0003191991199752, 'HP:0045027': 2.231477614568821e-03, 'Pathogenic': 0.454434532586802}

NAMES = {'HP:0033127': 'Musculoskeletal',
 'HP:0040064': 'Limbs',
 'HP:0000707': 'Nervous',
 'HP:0001939': 'Metabolism/homeostasis',
 'HP:0000152': 'Head/neck',
 'HP:0001626': 'Cardiovascular',
 'HP:0000119': 'Genitourinary',
 'HP:0000478': 'Eye',
 'HP:0002715': 'Immune',
 'HP:0001574': 'Integument',
 'HP:0001871': 'Blood/blood-forming tissues',
 'HP:0025031': 'Digestive',
 'HP:0002664': 'Neoplasm',
 'HP:0002086': 'Respiratory',
 'HP:0000818': 'Endocrine',
 'HP:0000598': 'Ear',
 'HP:0025354': 'Cellular',
 'HP:0001197': 'Prenatal development/birth',
 'HP:0001507': 'Growth',
 'HP:0025142': 'Constitutional',
 'HP:0000769': 'Breast',
 'HP:0001608': 'Voice',
 'HP:0045027': 'Thoracic cavity',
 'Pathogenic': 'Pathogenic'}

# Load additional annotations and training data
X_train_cols = pd.read_csv(PROJECT_DIR + '/data/train_columns.csv', header=None)[0].tolist()
genedata = pd.read_parquet(PROJECT_DIR + '/data/gene_data.pq')

def generate_train(X_train_columns, predictor):
    selected_features = joblib.load(PROJECT_DIR + '/selected_features/' + PREFIX + predictor + '_selected_features.joblib')
    columns = [c for c in X_train_columns if c != 'SYMBOL'] + selected_features
    return columns

X_train_br = generate_train(X_train_cols, 'br')
X_train_brsampling = generate_train(X_train_cols, 'brsampling')
X_train_rakel = generate_train(X_train_cols, 'rakel')
X_train_rakelsampling = generate_train(X_train_cols, 'rakelsampling')
X_train_lp = generate_train(X_train_cols, 'lp')
X_train_lpsampling = generate_train(X_train_cols, 'lpsampling')

paths = {'br': PROJECT_DIR + '/models/' + PREFIX + 'binaryrelevance.joblib',
         'brsampling': PROJECT_DIR + '/models/' + PREFIX + 'binaryrelevance_sampling.joblib',
         'rakel': PROJECT_DIR + '/models/' + PREFIX + 'rakeld.joblib',
         'rakelsampling': PROJECT_DIR + '/models/' + PREFIX + 'rakeld_sampling.joblib',
         'lp': PROJECT_DIR + '/models/' + PREFIX + 'labelpowerset.joblib',
         'lpsampling': PROJECT_DIR + '/models/' + PREFIX + 'labelpowerset_sampling.joblib'}

def predict(payload):
    predictor, df = payload
    df = df.replace('', None)
    df = df.replace('-', None)

    selected_features = joblib.load(PROJECT_DIR + '/selected_features/' + PREFIX + predictor + '_selected_features.joblib')
    preprocessor = joblib.load(PROJECT_DIR + '/preprocessors/' + PREFIX + predictor + '_preprocessor.joblib')

    df_temp = df.merge(genedata[selected_features], on='SYMBOL', how='left')
    if predictor == 'br':
        df_temp = df_temp[X_train_br]
    elif predictor == 'brsampling':
        df_temp = df_temp[X_train_brsampling]
    elif predictor == 'rakel':
        df_temp = df_temp[X_train_rakel]
    elif predictor == 'rakelsampling':
        df_temp = df_temp[X_train_rakelsampling]
    elif predictor == 'lp':
        df_temp = df_temp[X_train_lp]
    elif predictor == 'lpsampling':
        df_temp = df_temp[X_train_lpsampling]

    proc_data = preprocessor.transform(df_temp)
    model = joblib.load(paths[predictor])
    try:
        model.set_params(**{'n_jobs': 1})
    except:
        pass
    if predictor in ['br', 'brsampling']:
        model.estimator.set_params(**{'n_jobs': 1})
    elif predictor in ['rakel', 'rakelsampling']:
        model.base_classifier.set_params(**{'n_jobs': 1})
    else:
        model.classifier.set_params(**{'n_jobs': 1})

    pred = model.predict_proba(proc_data)

    if predictor in ['br', 'brsampling']:
        pred = np.array([yp[:, 1] for yp in pred]).T
    elif predictor in ['rakel', 'rakelsampling']:
        pred = pred.toarray()
    return pred

infiles = sys.argv[1:]
outnames = [f.replace('.pq', '_preds.csv') for f in infiles]

for infile, outname in zip(infiles, outnames):
    try:
        parquet_file = pq.ParquetFile(infile)
    except:
        continue

    for batch in tqdm(parquet_file.iter_batches(batch_size=100000)):
        df = batch.to_pandas()
        df = df.rename(columns={'Grantham_x':'Grantham', 'bStatistic_x': 'bStatistic'})
        df = df.reset_index(drop=True)

        with mp.pool.Pool(6) as p:
            all_pred = p.map(predict, zip(['br', 'brsampling', 'rakel', 'rakelsampling', 'lp', 'lpsampling'], itertools.repeat(df, 6)))

        out = pd.DataFrame(np.mean(all_pred, axis=0), columns=OUT_COLUMNS)
        crisp = []
        for out_i, (_, r) in enumerate(out.iterrows()):
            rc = [NAMES[l] for p, l in zip(r, out.columns) if p >= CUTOFFS[l]]
            if not rc:
                rc = ['Benign']
            crisp.append(','.join(rc))
        out['V2P_predicted_phenotypes'] = crisp
        out['ID'] = df['ID']
        out.columns = [NAMES[c] if c in NAMES.keys() else c for c in out.columns]
        if os.path.exists(outname):
            out.to_csv(outname, index=None, mode='a', header=None)
        else:
            out.to_csv(outname, index=None)
