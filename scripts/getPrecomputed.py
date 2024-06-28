import sqlite3
import genomicsqlite
import sys
import pandas as pd
import tabix
import os

PROJECT_DIR = os.environ['V2P_DIR']
DATABASE_LOCATION = sys.argv[3]
NUCS = ['A','T','G','C', 'a', 't', 'g', 'c']

TRAIN_CUTOFFS = {'HP:0033127': 0.18985256120541105, 'HP:0040064': 0.1992395353086459, 'HP:0000707': 0.27763873171701137, 'HP:0001939': 0.2551450070191609, 'HP:0000152': 0.15904484246054185, 'HP:0001626': 0.23145503270126955, 'HP:0000119': 0.23078016638280274, 'HP:0000478': 0.24262252332577205, 'HP:0002715': 0.20793127139198309, 'HP:0001574': 0.2532074285138736, 'HP:0001871': 0.2516688747291588, 'HP:0025031': 0.2233320540843347, 'HP:0002664': 0.24690386081384114, 'HP:0002086': 0.15949115580936354, 'HP:0000818': 0.24717306234251663, 'HP:0000598': 0.22988786600039773, 'HP:0025354': 0.21765933751203598, 'HP:0001197': 0.0143026400691627, 'HP:0001507': 0.1909110691851463, 'HP:0025142': 0.044549904297502496, 'HP:0000769': 0.2063466967995205, 'HP:0001608': 0.00392180555716705, 'HP:0045027': 6.146926123748867e-03, 'Pathogenic': 0.454434532586802}
TEST_CUTOFFS = {'HP:0033127': 0.09636672984829231, 'HP:0040064': 0.02934571321747775, 'HP:0000707': 0.19613093774473234, 'HP:0001939': 0.0223904127687745, 'HP:0000152': 0.09796504479742005, 'HP:0001626': 0.02974164625117005, 'HP:0000119': 0.06618039735346559, 'HP:0000478': 0.0677496075849934, 'HP:0002715': 0.095849938545434, 'HP:0001574': 0.04389976814158695, 'HP:0001871': 0.07640677609666524, 'HP:0025031': 0.05631736513563655, 'HP:0002664': 0.2433363216373498, 'HP:0002086': 0.2087418508756681, 'HP:0000818': 0.04985664742061175, 'HP:0000598': 0.06204993227308815, 'HP:0025354': 0.07053955236430914, 'HP:0001197': 0.0030437757064105, 'HP:0001507': 0.0224615051575628, 'HP:0025142': 0.0144888249770768, 'HP:0000769': 0.08022135225114585, 'HP:0001608': 0.00033938205675165, 'HP:0045027': 2.3046411555773557e-03, 'Pathogenic': 0.4011821405315217}

NAMES = {'musculoskeletal': 'Musculoskeletal',
 'limbs': 'Limbs',
 'nervous': 'Nervous',
 'metabolism': 'Metabolism/homeostasis',
 'head': 'Head/neck',
 'cardiovascular': 'Cardiovascular',
 'genitourinary': 'Genitourinary',
 'eye': 'Eye',
 'immune': 'Immune',
 'integument': 'Integument',
 'blood': 'Blood/blood-forming tissues',
 'digestive': 'Digestive',
 'neoplasm': 'Neoplasm',
 'respiratory': 'Respiratory',
 'endocrine': 'Endocrine',
 'ear': 'Ear',
 'cellular': 'Cellular',
 'prenatal': 'Prenatal development/birth',
 'growth': 'Growth',
 'constitutional': 'Constitutional',
 'breast': 'Breast',
 'voice': 'Voice',
 'thoracic': 'Thoracic cavity',
 'Pathogenic': 'Pathogenic'}

NAMES_HPO = {'Musculoskeletal': 'HP:0033127',
 'Limbs': 'HP:0040064',
 'Nervous': 'HP:0000707',
 'Metabolism/homeostasis': 'HP:0001939',
 'Head/neck': 'HP:0000152',
 'Cardiovascular': 'HP:0001626',
 'Genitourinary': 'HP:0000119',
 'Eye': 'HP:0000478',
 'Immune': 'HP:0002715',
 'Integument': 'HP:0001574',
 'Blood/blood-forming tissues': 'HP:0001871',
 'Digestive': 'HP:0025031',
 'Neoplasm': 'HP:0002664',
 'Respiratory': 'HP:0002086',
 'Endocrine': 'HP:0000818',
 'Ear': 'HP:0000598',
 'Cellular': 'HP:0025354',
 'Prenatal development/birth': 'HP:0001197',
 'Growth': 'HP:0001507',
 'Constitutional': 'HP:0025142',
 'Breast': 'HP:0000769',
 'Voice': 'HP:0001608',
 'Thoracic cavity': 'HP:0045027',
 'Pathogenic': 'Pathogenic'}

ORDER = ['musculoskeletal', 'limbs', 'nervous', 'metabolism', 'head', 'cardiovascular', 'genitourinary', 'eye',
         'immune', 'integument', 'blood', 'digestive', 'neoplasm', 'respiratory', 'endocrine', 'ear', 'cellular',
         'prenatal', 'growth', 'constitutional', 'breast', 'voice', 'thoracic', 'Pathogenic', 'uid']

train_symbols = set(pd.read_csv(PROJECT_DIR + '/data/train_symbols.csv')['SYMBOL'].tolist())
tb = tabix.open(PROJECT_DIR + '/data/symbols.csv.gz')

data = pd.read_csv(sys.argv[1], sep='\t')
chroms = list(data['#CHROM'].unique())

tmpfile = sys.argv[2]

def get_crisp(df, symbols):
    crisp = []
    for index, (_, r) in enumerate(df.iterrows()):
        CUTOFFS = TRAIN_CUTOFFS if symbols[index] in train_symbols else TEST_CUTOFFS
        rc = [NAMES[l] for p, l in zip(r, df.columns) if l != 'uid' and p >= CUTOFFS[NAMES_HPO[NAMES[l]]]]
        if not rc:
            rc = ['Benign']
        crisp.append(','.join(rc))
    return crisp

def get_symbols(df):
    chrom = df['uid'].apply(lambda v: v.split('_')[0])
    pos = df['uid'].apply(lambda v: v.split('_')[1])
    symbols = []
    for c, p in zip(chrom, pos):
        q = c + ':' + p + '-' + p
        res = list(tb.querys(q))
        if not res:
            symbols.append(None)
        else:
            symbols.append(res[0][0])
    return symbols

def get_predictions(conn):
    cursor = conn.cursor()
    query = f'''
    SELECT {','.join(ORDER)} 
    FROM predictions 
    WHERE uid IN {repr(tuple(map(str, snps['ID'].tolist())))}
    '''
    cursor.execute(query)
    return cursor.fetchall()

for chrom in chroms:

    chrom_df = data.loc[data['#CHROM'] == chrom]
    snps = chrom_df.loc[((chrom_df['REF'].isin(NUCS)) & (chrom_df['ALT'].isin(NUCS)))]
    indels = chrom_df.loc[~chrom_df['ID'].isin(snps['ID'])]

    snp_dbconn = genomicsqlite.connect(
        DATABASE_LOCATION + 'chrom' + str(chrom) + '_compressed.db',
        read_only=True,
    )

    indel_dbconn = genomicsqlite.connect(
        DATABASE_LOCATION + 'gnomad_chrom' + str(chrom) + '_compressed.db',
        read_only=True,
    )

    snp_rows = pd.DataFrame(get_predictions(snp_dbconn), columns=ORDER)
    indel_rows = pd.DataFrame(get_predictions(indel_dbconn), columns=ORDER)

    snp_symbols = get_symbols(snp_rows)
    indel_symbols = get_symbols(indel_rows)
    
    snp_rows['V2P_predicted_phenotypes'] = get_crisp(snp_rows, snp_symbols)
    indel_rows['V2P_predicted_phenotypes'] = get_crisp(indel_rows, indel_symbols)

    snp_rows.to_csv(tmpfile, index=None, header=None, mode='a+')
    indel_rows.to_csv(tmpfile, index=None, header=None, mode='a+')

found = pd.read_csv(tmpfile, header=None)
data = data.loc[~data['ID'].isin(found[found.columns[-1]])]
data.to_csv(sys.argv[1], index=None, sep='\t')
