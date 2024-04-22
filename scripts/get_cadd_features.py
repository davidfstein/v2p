import sys
import pandas as pd
import numpy as np
import sqlite3
import genomicsqlite
from tqdm import tqdm
from column_map import cadd_feat_map

CADD_COLS = ['Chrom','Pos','Ref','Alt','ConsScore','GC','CpG','motifECount','motifEName','motifEHIPos','motifEScoreChng','GeneID','Dst2Splice','Dst2SplType','minDistTSS','minDistTSE','priPhCons','mamPhCons','verPhCons','priPhyloP','mamPhyloP','verPhyloP','bStatistic','targetScan','mirSVR_Score','mirSVR_E','mirSVR_Aln','cHmm_E1','cHmm_E2','cHmm_E3','cHmm_E4','cHmm_E5','cHmm_E6','cHmm_E7','cHmm_E8','cHmm_E9','cHmm_E10','cHmm_E11','cHmm_E12','cHmm_E13','cHmm_E14','cHmm_E15','cHmm_E16','cHmm_E17','cHmm_E18','cHmm_E19','cHmm_E20','cHmm_E21','cHmm_E22','cHmm_E23','cHmm_E24','cHmm_E25','GerpRS','GerpRSpval','GerpN','GerpS','tOverlapMotifs','motifDist','EncodeH3K4me1_sum','EncodeH3K4me1_max','EncodeH3K4me2_sum','EncodeH3K4me2_max','EncodeH3K4me3_sum','EncodeH3K4me3_max','EncodeH3K9ac_sum','EncodeH3K9ac_max','EncodeH3K9me3_sum','EncodeH3K9me3_max','EncodeH3K27ac_sum','EncodeH3K27ac_max','EncodeH3K27me3_sum','EncodeH3K27me3_max','EncodeH3K36me3_sum','EncodeH3K36me3_max','EncodeH3K79me2_sum','EncodeH3K79me2_max','EncodeH4K20me1_sum','EncodeH4K20me1_max','EncodeH2AFZ_sum','EncodeH2AFZ_max','EncodeDNase_sum','EncodeDNase_max','EncodetotalRNA_sum','EncodetotalRNA_max','Grantham','SpliceAI_acc_gain','SpliceAI_acc_loss','SpliceAI_don_gain','SpliceAI_don_loss','MMSp_acceptorIntron','MMSp_acceptor','MMSp_exon','MMSp_donor','MMSp_donorIntron','Dist2Mutation','Freq100bp','Rare100bp','Sngl100bp','Freq1000bp','Rare1000bp','Sngl1000bp','Freq10000bp','Rare10000bp','Sngl10000bp','EnsembleRegulatoryFeature','dbscSNV_ada_score','dbscSNV_rf_score','RemapOverlapTF','RemapOverlapCL']
SELECT_QUERY = 'SELECT ' + ','.join(CADD_COLS) + ' FROM cadd WHERE Chrom="%s" and Pos="%s" and Ref="%s" and Alt="%s" and GeneID="%s"'
BACKUP_SELECT_QUERY = 'SELECT ' + ','.join(CADD_COLS) + ' FROM cadd WHERE Chrom="%s" and Pos="%s" and Ref="%s" and Alt="%s"'

def backup_query(r, cursor):
    _ = cursor.execute(BACKUP_SELECT_QUERY%(r['CHROM'], r['POS'], r['REF'], r['ALT']))
    rows = cursor.fetchall()
    if not len(rows):
        return []
    if len(rows) > 1:
        cs = np.array([sum(v is not None for v in row) for row in rows])
        return [rows[np.argmax(cs)]]
    else:
        return rows

def get_cadd_features(path, threads=1):
    data = pd.read_parquet(path, columns=['ID', 'CHROM', 'POS', 'REF', 'ALT', 'Gene'])
    data.index = data['CHROM'].astype(str)
    chroms = [str(i+1) for i in range(22)]
    chroms.extend(['Y', 'X'])
    out = []
    for chrom in chroms:
        conn = genomicsqlite.connect(
                '/cadddb/cadd' + str(chrom) + '_compressed.db',
                read_only=True,
                **{'isolation_level': 'DEFERRED', 'threads': int(threads), 'immutable': True}
        )
        cursor = conn.cursor()
        try:
            # Optimizations
            cursor.execute('pragma journal_mode = WAL;')
            cursor.execute('pragma synchronous = normal;')
            cursor.execute('pragma temp_store = memory;')
            cursor.execute('pragma mmap_size = 30000000000;')
            if chrom not in data.index:
                continue
            df = data.loc[[chrom]]
            for _, r in tqdm(df.iterrows(), total=df.shape[0]):
                _ = cursor.execute(SELECT_QUERY%(r['CHROM'], r['POS'], r['REF'], r['ALT'], r['Gene']))
                rows = cursor.fetchall()
                if not len(rows):
                    rows = backup_query(r, cursor)
                    if not len(rows):
                        continue
                if len(rows) != 1:
                    raise ValueError(rows, r)
                out.append([r['ID']] + list(rows[0]))
        except Exception as e:
            cursor.close()
            raise e
        finally:
            cursor.close()
    return pd.DataFrame(out, columns=['ID'] + [cadd_feat_map[c] if c in cadd_feat_map.keys() else c for c in CADD_COLS])

if __name__ == "__main__":
    out = get_cadd_features(sys.argv[1], sys.argv[2])
