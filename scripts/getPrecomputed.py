import sqlite3
import genomicsqlite
import sys
import pandas as pd
import os

PROJECT_DIR = os.environ['V2P_DIR']
DATABASE_LOCATION = sys.argv[3]
NUCS = ['A','T','G','C', 'a', 't', 'g', 'c']

CUTOFFS = {'HP:0033127': 0.13375073118275554, 'HP:0040064': 0.1027121652336587, 'HP:0000707': 0.1280226940316797, 'HP:0001939': 0.13289699594521331, 'HP:0000152': 0.10556380274414325, 'HP:0001626': 0.1553588637596261, 'HP:0000119': 0.1332160134368378, 'HP:0000478': 0.20826082412161753, 'HP:0002715': 0.15873701583024458, 'HP:0001574': 0.16293394524364474, 'HP:0001871': 0.13980646528040241, 'HP:0025031': 0.1296962674378072, 'HP:0002664': 0.10575180529819386, 'HP:0002086': 0.02704949652042875, 'HP:0000818': 0.12326556052836206, 'HP:0000598': 0.12354089044413866, 'HP:0025354': 0.23111810436978011, 'HP:0001197': 0.0020001575650973, 'HP:0001507': 0.21537008781720823, 'HP:0025142': 0.0988696186665396, 'HP:0000769': 0.13120347198135124, 'HP:0001608': 0.0003191991199752, 'HP:0045027': 2.231477614568821e-03, 'Pathogenic': 0.454434532586802}

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
 'Pathogenic': 'Pathogenic',
 'uid': 'ID'}

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
OUT_ORDER = ['ID','V2P_predicted_phenotypes','Musculoskeletal','Limbs','Nervous','Metabolism/homeostasis','Head/neck','Cardiovascular','Genitourinary',
             'Eye','Immune','Integument','Blood/blood-forming tissues','Digestive','Neoplasm','Respiratory','Endocrine',
             'Ear','Cellular','Prenatal development/birth','Growth','Constitutional','Breast','Voice','Thoracic cavity',
             'Pathogenic']

data = pd.read_csv(sys.argv[1], sep='\t', low_memory=False)
data['#CHROM'] = data['#CHROM'].apply(lambda v: str(v).replace('chr',''))
data['ID'] = data.apply(lambda r: '_'.join([str(r['#CHROM']), str(int(r['POS'])),r['REF'], r['ALT']]), axis=1) 
chroms = list(data['#CHROM'].unique())

tmpfile = sys.argv[2]

def get_crisp(df):
    crisp = []
    for _, r in df.iterrows():
        rc = [l for p, l in zip(r, df.columns) if l != 'ID' and p >= CUTOFFS[NAMES_HPO[l]]]
        if not rc:
            rc = ['Benign']
        crisp.append(','.join(rc))
    return crisp

def get_predictions(conn, variants):
    cursor = conn.cursor()
    query = f'''
    SELECT {','.join(ORDER)} 
    FROM predictions 
    WHERE uid IN {repr(tuple(map(str, variants['ID'].tolist() + [''])))}
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

    snp_rows = pd.DataFrame(columns=[NAMES[c] for c in ORDER])
    indel_rows = pd.DataFrame(columns=[NAMES[c] for c in ORDER])
    if not snps.empty:
        snp_rows = pd.DataFrame(get_predictions(snp_dbconn, snps), columns=[NAMES[c] for c in ORDER])
    if not indels.empty:
        indel_rows = pd.DataFrame(get_predictions(indel_dbconn, indels), columns=[NAMES[c] for c in ORDER])

    snp_rows['V2P_predicted_phenotypes'] = get_crisp(snp_rows)
    indel_rows['V2P_predicted_phenotypes'] = get_crisp(indel_rows)

    file_populated = os.path.exists(tmpfile) and os.stat(tmpfile).st_size > 0
    header = None if file_populated else True
    mode = 'a+' if file_populated else 'w+'
    snp_rows[OUT_ORDER].to_csv(tmpfile, index=None, header=header, mode=mode)
    indel_rows[OUT_ORDER].to_csv(tmpfile, index=None, header=None, mode='a+')

found = pd.read_csv(tmpfile, usecols=['ID'])
data = data.loc[~data['ID'].isin(found['ID'])]
outname = sys.argv[1].split('/')[-1].replace('.vcf', '_novel.vcf')
data.to_csv(outname, index=None, sep='\t')
