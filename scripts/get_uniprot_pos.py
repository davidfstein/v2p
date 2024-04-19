import pandas as pd
from tqdm import tqdm
import sys
import sqlite3

def main():
    data = pd.read_parquet(sys.argv[1], columns=['POS'])
    ids = pd.read_csv('lookup_ids.csv')
    data['proteinID'] = ids['proteinID']
    select_query = 'SELECT aa_pos FROM uniprot_genomic_position WHERE proteinID=%s and genomic_pos=%s'
    out = []
    for _, r in tqdm(data.iterrows(), total=data.shape[0]):
        if pd.isna(r['proteinID']):
            out.append(None)
            continue
        
        qry = select_query%(r['proteinID'], r['POS'])
        _ = cursor.execute(qry)
        rows = cursor.fetchall()
        if len(rows) == 0:
            out.append(None)
            continue
        if len(rows) != 1:
            raise ValueError(rows, r)
        out.append(rows[0][0])
    ids['uniprot_pos'] = out
    return ids

if __name__ == '__main__':
    conn = sqlite3.connect('/home/hpo.db', isolation_level='DEFERRED')
    cursor = conn.cursor()
    try:
        # Optimizations
        cursor.execute('pragma journal_mode = WAL;')
        cursor.execute('pragma synchronous = normal;')
        cursor.execute('pragma temp_store = memory;')
        cursor.execute('pragma mmap_size = 30000000000;')
        ids = main()
        ids.to_csv('lookup_ids.csv', index=None)
    except Exception as e:
        raise e
    finally:
        cursor.close()