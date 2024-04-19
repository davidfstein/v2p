import pandas as pd
import sqlite3
from tqdm import tqdm
from column_map import misc_gene_map

GET_COLUMNS_QUERY="SELECT name FROM PRAGMA_TABLE_INFO('%s');"

def get_misc_gene_features(data, cursor):
    _ = cursor.execute(GET_COLUMNS_QUERY%"misc_gene_features")
    rows = cursor.fetchall()
    # skip first which is geneID
    columns = [row[0] for row in rows[1:]]
    select_query = 'SELECT ' + ','.join(columns) + ' FROM misc_gene_features WHERE geneID=%s'
    out = []
    for gid in tqdm(data['geneID']):
        if pd.isna(gid):
            out.append([None] * len(columns))
            continue
        
        qry = select_query%gid
        _ = cursor.execute(qry)
        rows = cursor.fetchall()
        if len(rows) == 0:
            out.append([None] * len(columns))
            continue
        if len(rows) != 1:
            raise ValueError(rows, gid)
        out.append(rows[0])
    return pd.DataFrame(out, columns=[misc_gene_map[c] if c in misc_gene_map.keys() else c for c in columns])

def get_geneID_features():
    conn = sqlite3.connect('/home/hpo.db', isolation_level='DEFERRED')
    cursor = conn.cursor()
    try:
        # Optimizations
        cursor.execute('pragma journal_mode = WAL;')
        cursor.execute('pragma synchronous = normal;')
        cursor.execute('pragma temp_store = memory;')
        cursor.execute('pragma mmap_size = 30000000000;')
        ids = pd.read_csv('lookup_ids.csv')
        return get_misc_gene_features(ids, cursor)
    except Exception as e:
        print(e)
    finally:
        cursor.close()

if __name__ == "__main__":
    get_geneID_features()