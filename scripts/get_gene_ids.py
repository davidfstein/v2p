import pandas as pd
from tqdm import tqdm
import sqlite3
import sys

def main():
    data = pd.read_parquet(sys.argv[1], columns=['SYMBOL', 'Gene'])
    data = data.replace('.', None)
    gene_symbol_query = 'SELECT geneID FROM gene WHERE SYMBOL="%s" and ensemblID="%s"'
    symbol_query = 'SELECT geneID FROM gene WHERE SYMBOL="%s"'
    gene_query = 'SELECT geneID FROM gene WHERE ensemblID="%s"'

    seen = {}
    qry=None
    gids = []
    for _, r in tqdm(data.iterrows(), total=data.shape[0]):
        identifier = None
        if not pd.isna(r['Gene']) and not pd.isna(r['SYMBOL']):
            qry = gene_symbol_query%(r['SYMBOL'], r['Gene'])
            identifier = "_".join([r['SYMBOL'], r['Gene']])
        elif not pd.isna(r['SYMBOL']):
            qry = symbol_query%r['SYMBOL']
            identifier = "_".join([r['SYMBOL']])
        elif not pd.isna(r['Gene']):
            qry = gene_query%r['Gene']
            identifier = "_".join([r['Gene']])
        else:
            qry = None
        if r['Gene'] == '.' and r['SYMBOL'] == '.':
            qry = None
        if not qry:
            gids.append(None)
            continue

        if identifier in seen.keys():
            gids.append(seen[identifier])
            continue
    
        _ = cursor.execute(qry)
        rows = cursor.fetchall()
        if len(rows) != 1 or len(rows[0]) != 1:
            print(rows, '\n', r)
            raise ValueError()
        seen[identifier] = rows[0][0]
        gids.append(rows[0][0])
    pd.DataFrame(gids,columns=['geneID']).to_csv(sys.argv[2], index=None) 

if __name__ == '__main__':
    conn = sqlite3.connect('/home/hpo.db', isolation_level='DEFERRED')
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
