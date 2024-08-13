import pandas as pd
from tqdm import tqdm
import sqlite3
import sys

def main():
    data = pd.read_parquet(sys.argv[1], columns=['ENSP', 'SWISSPROT'])
    data = data.replace('.', None)
    data['UNIPROT'] = data['SWISSPROT'].apply(lambda v: v.split('.')[0] if not pd.isna(v) else v)
    gids = pd.read_csv('lookup_ids.csv')
    ensp_uniprot_query = 'SELECT proteinID FROM protein WHERE geneID=%s and uniprot="%s" and ensp="%s"'
    uniprot_query = 'SELECT proteinID FROM protein WHERE geneID=%s and uniprot="%s"'
    ensp_query = 'SELECT proteinID FROM protein WHERE geneID=%s and ensp="%s"'

    seen = {}
    qry=None
    pids = []
    for i, r in tqdm(data.iterrows(), total=data.shape[0]):
        identifier = None
        if not pd.isna(r['UNIPROT']) and not pd.isna(r['ENSP']):
            qry = ensp_uniprot_query%(gids['geneID'][i], r['UNIPROT'], r['ENSP'])
            identifier = "_".join([str(gids['geneID'][i]), r['UNIPROT'], r['ENSP']])
        elif not pd.isna(r['UNIPROT']):
            qry = uniprot_query%(gids['geneID'][i], r['UNIPROT'])
            identifier = "_".join([str(gids['geneID'][i]), r['UNIPROT']])
        elif not pd.isna(r['ENSP']):
            qry = ensp_query%(gids['geneID'][i], r['ENSP'])
            identifier = "_".join([str(gids['geneID'][i]), r['ENSP']])
        else:
            qry = None
        if r['UNIPROT'] in ['.',''] and r['ENSP'] == '.':
            qry = None
        if not qry:
            pids.append(None)
            continue

        if identifier in seen.keys():
            pids.append(seen[identifier])
            continue
    
        _ = cursor.execute(qry)
        rows = cursor.fetchall()
        if len(rows) != 1 or len(rows[0]) != 1:
            print(identifier, rows, r, sep='\n')
            raise ValueError()
        seen[identifier] = rows[0][0]
        pids.append(rows[0][0])
    gids['proteinID'] = pids
    gids.to_csv('lookup_ids.csv', index=None)

if __name__ == '__main__':
    conn = sqlite3.connect('/home/myuser/work/hpo.db', isolation_level='DEFERRED')
    cursor = conn.cursor()
    try:
        # Optimizations
        cursor.execute('pragma journal_mode = WAL;')
        cursor.execute('pragma synchronous = normal;')
        cursor.execute('pragma temp_store = memory;')
        cursor.execute('pragma mmap_size = 30000000000;')
        main()
    finally:
        cursor.close()
