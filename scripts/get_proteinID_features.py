import pandas as pd
import sqlite3
import sys
from multiprocessing import Process, Manager
from tqdm import tqdm
from column_map import uniprot_aa_map, uniprot_mut_features_map

GET_COLUMNS_QUERY="SELECT name FROM PRAGMA_TABLE_INFO('%s');"

# names are the same
def get_misc_protein_features(data, cursor, result_dict=None):
    _ = cursor.execute(GET_COLUMNS_QUERY%"misc_protein_features")
    rows = cursor.fetchall()
    # skip first which is proteinID
    columns = [row[0] for row in rows[1:]]
    select_query = 'SELECT ' + ','.join(columns) + ' FROM misc_protein_features WHERE proteinID=%s'
    out = []
    for pid in tqdm(data['proteinID']):
        if pd.isna(pid):
            out.append([None] * len(columns))
            continue
        
        qry = select_query%pid
        _ = cursor.execute(qry)
        rows = cursor.fetchall()
        if len(rows) == 0:
            out.append([None] * len(columns))
            continue
        if len(rows) != 1:
            tr = rows[0]
            if all([r == tr for r in rows]):
                rows = [rows[0]]
            else:
                raise ValueError(rows, pid)
        out.append(rows[0])
    out = pd.DataFrame(out, columns=columns)
    if result_dict is not None:
        result_dict['mpf'] = out
    return out 

# names are the same
def get_ensp_aa_feat(data, cursor, result_dict=None):
    _ = cursor.execute(GET_COLUMNS_QUERY%"ensp_aa_feat")
    rows = cursor.fetchall()
    # skip first three which are proteinID, aa_pos, and aa_name
    columns = [row[0] for row in rows[3:]]
    select_query = 'SELECT ' + ','.join(columns) + ' FROM ensp_aa_feat WHERE proteinID=%s and aa_pos=%s'
    out = []
    for _, r in tqdm(data.iterrows(), total=data.shape[0]):
        if pd.isna(r['proteinID']) or pd.isna(r['AA_pos']):
            out.append([None] * len(columns))
            continue

        qry = select_query%(r['proteinID'], r['AA_pos'])
        _ = cursor.execute(qry)
        rows = cursor.fetchall()
        if len(rows) == 0:
            out.append([None] * len(columns))
            continue
        if len(rows) != 1:
            tr = rows[0]
            if all([r == tr for r in rows]):
                rows = [rows[0]]
            else:
                raise ValueError(rows, r)
        out.append(rows[0])
    out = pd.DataFrame(out, columns=columns)
    if result_dict  is not None:
        result_dict['eaf'] = out
    return out 

def get_ensp_mut_features(data, cursor, result_dict=None):
    select_query = 'SELECT S_DDG_SEQ FROM ensp_mut_features WHERE proteinID=%s and aa_pos=%s and alt_aa_name="%s"'
    out = []
    for _, r in tqdm(data.iterrows(), total=data.shape[0]):
        if pd.isna(r['proteinID']) or pd.isna(r['AA_pos']) or pd.isna(r['Alt_aa']):
            out.append([None])
            continue

        qry = select_query%(r['proteinID'], r['AA_pos'], r['Alt_aa'])
        _ = cursor.execute(qry)
        rows = cursor.fetchall()
        if len(rows) == 0:
            out.append([None])
            continue
        if len(rows) != 1:
            raise ValueError(rows, r)
        out.append(rows[0])
    out = pd.DataFrame(out, columns=['S_DDG[SEQ]'])
    if result_dict is not None:
        result_dict['emf'] = out
    return out 

def get_uniprot_aa_features(data, cursor, result_dict=None):
    _ = cursor.execute(GET_COLUMNS_QUERY%"uniprot_aa_feat")
    rows = cursor.fetchall()
    # skip first three which are proteinID, aa_pos, and genomic_pos
    columns = [row[0] for row in rows[3:]]
    select_query = 'SELECT ' + ','.join(columns) + ' FROM uniprot_aa_feat WHERE proteinID=%s and aa_pos=%s'
    out = []
    for _, r in tqdm(data.iterrows(), total=data.shape[0]):
        if pd.isna(r['proteinID']) or pd.isna(r['uniprot_pos']):
            out.append([None] * len(columns))
            continue

        qry = select_query%(r['proteinID'], r['uniprot_pos'])
        _ = cursor.execute(qry)
        rows = cursor.fetchall()
        if len(rows) == 0:
            out.append([None] * len(columns))
            continue
        if len(rows) != 1:
            raise ValueError(rows, r)
        out.append(rows[0])
    out = pd.DataFrame(out, columns=[uniprot_aa_map[c] if c in uniprot_aa_map.keys() else c for c in columns])
    if result_dict is not None:
        result_dict['uaf'] = out
    return out 

def get_uniprot_mut_features(data, cursor, result_dict=None):
    _ = cursor.execute(GET_COLUMNS_QUERY%"uniprot_mut_features")
    rows = cursor.fetchall()
    # skip first four which are proteinID, aa_pos, ref_aa_name, and alt_aa_name
    columns = [row[0] for row in rows[4:]]
    select_query = 'SELECT ' + ','.join(columns) + ' FROM uniprot_mut_features WHERE proteinID=%s and aa_pos=%s and alt_aa_name="%s"'
    out = []
    for _, r in tqdm(data.iterrows(), total=data.shape[0]):
        if pd.isna(r['proteinID']) or pd.isna(r['uniprot_pos']) or pd.isna(r['Alt_aa']):
            out.append([None] * len(columns))
            continue

        qry = select_query%(r['proteinID'], r['uniprot_pos'], r['Alt_aa'])
        _ = cursor.execute(qry)
        rows = cursor.fetchall()
        if len(rows) == 0:
            out.append([None] * len(columns))
            continue
        if len(rows) != 1:
            raise ValueError(rows, r)
        out.append(rows[0])
    out = pd.DataFrame(out, columns=[uniprot_mut_features_map[c] if c in uniprot_mut_features_map.keys() else c for c in columns])
    if result_dict is not None:
        result_dict['umf'] = out
    return out 

def get_proteinID_features(path, n_cores=1):
    conn = sqlite3.connect('/home/hpo.db', isolation_level='DEFERRED')
    cursor = conn.cursor()
    try:
        # Optimizations
        cursor.execute('pragma journal_mode = WAL;')
        cursor.execute('pragma synchronous = normal;')
        cursor.execute('pragma temp_store = memory;')
        cursor.execute('pragma mmap_size = 30000000000;')
        ids = pd.read_csv('lookup_ids.csv')
        data = pd.read_parquet(path, columns=['AA_pos', 'Alt_aa'])
        data['proteinID'] = ids['proteinID']
        data['uniprot_pos'] = ids['uniprot_pos']
        f1 = []
        f2 = []
        f3 = []
        f4 = []
        f5 = []
        if int(n_cores) >= 5:
            manager = Manager()
            result_dict = manager.dict()
            mpf_proc = Process(target=get_misc_protein_features, args=(data, cursor), kwargs={'result_dict': result_dict})
            eaf_proc = Process(target=get_ensp_aa_feat, args=(data, cursor), kwargs={'result_dict': result_dict})
            emf_proc = Process(target=get_ensp_mut_features, args=(data, cursor), kwargs={'result_dict': result_dict})
            uaf_proc = Process(target=get_uniprot_aa_features, args=(data, cursor), kwargs={'result_dict': result_dict})
            umf_proc = Process(target=get_uniprot_mut_features, args=(data, cursor), kwargs={'result_dict': result_dict})
            procs = [mpf_proc, eaf_proc, emf_proc, uaf_proc, umf_proc]
            for proc in procs:
                proc.start()
            for proc in procs:
                proc.join()
            f1 = result_dict['mpf']
            f2 = result_dict['eaf']
            f3 = result_dict['emf']
            f4 = result_dict['uaf']
            f5 = result_dict['umf']
        else:
            f1 = get_misc_protein_features(data, cursor)
            f2 = get_ensp_aa_feat(data, cursor)
            f3 = get_ensp_mut_features(data, cursor)
            f4 = get_uniprot_aa_features(data, cursor)
            f5 = get_uniprot_mut_features(data, cursor)
        assert(f1.shape[0] == f2.shape[0] == f3.shape[0] == f4.shape[0] == f5.shape[0])
        return pd.concat([f1, f2, f3, f4, f5], axis=1)
    finally:
        cursor.close()

if __name__ == "__main__":
    get_proteinID_features(sys.argv[1], int(sys.argv[2]))
