import pandas as pd
import sqlite3
from tqdm import tqdm
import sys

GET_COLUMNS_QUERY="SELECT name FROM PRAGMA_TABLE_INFO('%s');"

def get_snvbox_aa_features(data, cursor):
    _ = cursor.execute(GET_COLUMNS_QUERY%"snvbox_aa_features")
    rows = cursor.fetchall()
    #skip first two which are WildType and Mut 
    columns = [row[0] for row in rows[2:]]
    select_query = 'SELECT ' + ','.join(columns) + ' FROM snvbox_aa_features WHERE WildType="%s" and Mut="%s"'
    out = []
    for _, r in tqdm(data.iterrows(), total=data.shape[0]):
        if pd.isna(r['Ref_aa']) or pd.isna(r['Alt_aa']):
            out.append([None] * len(columns))
            continue
        
        qry = select_query%(r['Ref_aa'], r['Alt_aa'])
        _ = cursor.execute(qry)
        rows = cursor.fetchall()
        if len(rows) == 0:
            out.append([None] * len(columns))
            continue
        if len(rows) != 1:
            raise ValueError(rows, r)
        out.append(rows[0])
    return pd.DataFrame(out, columns=columns)

def get_aa_features(path):
    conn = sqlite3.connect('/home/myuser/work/hpo.db', isolation_level='DEFERRED')
    cursor = conn.cursor()
    try:
        # Optimizations
        cursor.execute('pragma journal_mode = WAL;')
        cursor.execute('pragma synchronous = normal;')
        cursor.execute('pragma temp_store = memory;')
        cursor.execute('pragma mmap_size = 30000000000;')
        data = pd.read_parquet(path, columns=['Ref_aa', 'Alt_aa'])
        return get_snvbox_aa_features(data, cursor)
    finally:
        cursor.close()

if __name__ == "__main__":
    get_aa_features(sys.argv[1])
