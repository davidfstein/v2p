import sys
import math
import multiprocessing
import psutil
import resource
import time
import signal
import pandas as pd
import numpy as np
import scipy.io as sio
from numpy.linalg import inv
from numpy import linalg as LA
from sklearn.metrics.pairwise import cosine_similarity as cossim
from numpy import count_nonzero
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.compose import ColumnTransformer, make_column_selector
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import SelectFromModel, VarianceThreshold, SelectFdr, f_classif, chi2, mutual_info_classif
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.preprocessing import MinMaxScaler, OneHotEncoder, OrdinalEncoder
from sklearn.decomposition import PCA
from multiprocessing.pool import Pool
from multiprocessing import Manager
from functools import partial
import warnings
try:
    from imblearn.pipeline import Pipeline
    from imblearn.under_sampling import  RandomUnderSampler
    from imblearn.over_sampling import SMOTE, RandomOverSampler
except:
    from sklearn.pipeline import Pipeline
    print('Warning: imbalanced-learn not available', file=sys.stderr)
import random
from tqdm import tqdm
from lightgbm import LGBMClassifier

NEGONE_FEATURES = ['BLOSUM62','ProteinLengthChange','MaxEntScan_alt','MaxEntScan_diff','MaxEntScan_ref','ada_score','rf_score','rel_cDNA_pos','rel_CDS_pos','rel_prot_pos','Phosphorylation','Acetylation','Methylation','Ubiquitination','Glycosylation','PTM','RSA','ASA','AF_Relative_ASA','IUPRED2','ANCHOR2','A3D_SCORE','n_contacts','distance_com','concavity_score','S_DDG[SEQ]','S_DDG[3D]','hgmd_mutcount','gnomsingle_mutcount','gnom_mutcount','AF_confidence','isHomomultimer','num_interactions','ppi_combined_0','ppi_combined_1','ppi_combined_2','ppi_combined_3','ppi_combined_4','ppi_combined_5','ppi_combined_6','ppi_combined_7','ppi_combined_8','ppi_combined_9','ppi_combined_10','ppi_combined_11','ppi_combined_12','ppi_combined_13','ppi_combined_14','ppi_combined_15','ppi_combined_16','ppi_combined_17','ppi_combined_18','ppi_combined_19','ppi_combined_20','ppi_combined_21','ppi_combined_22','ppi_combined_23','ppi_combined_24','ppi_combined_25','ppi_combined_26','ppi_combined_27','ppi_combined_28','ppi_combined_29','ppi_combined_30','ppi_combined_31','ppi_combined_32','ppi_combined_33','ppi_combined_34','ppi_combined_35','ppi_combined_36','ppi_combined_37','ppi_combined_38','ppi_combined_39','ppi_combined_40','ppi_combined_41','ppi_combined_42','ppi_combined_43','ppi_combined_44','ppi_combined_45','ppi_combined_46','ppi_combined_47','ppi_combined_48','ppi_combined_49','ppi_combined_50','ppi_combined_51','ppi_combined_52','ppi_combined_53','ppi_combined_54','ppi_combined_55','ppi_combined_56','ppi_combined_57','ppi_combined_58','ppi_combined_59','ppi_combined_60','ppi_combined_61','ppi_combined_62','ppi_combined_63','DRNApredDNAscore_aa','ASAquick_normscore_aa','ASAquick_rawscore_aa','DFLpredScore_aa','DRNApredRNAscore_aa','DisoDNAscore_aa','DisoPROscore_aa','DisoRNAscore_aa','MMseq2_conservation_level_aa','MMseq2_conservation_score_aa','MoRFchibiScore_aa','PSIPRED_helix_aa','PSIPRED_strand_aa','SCRIBERscore_aa','SignalP_score_aa','PHOSPHORYLATION','ACETYLATION','UBIQUITINATION','S-NITROSYLATION','N-GLYCOSYLATION','METHYLATION','O-GLYCOSYLATION','MYRISTOYLATION','C-GLYCOSYLATION','SUMOYLATION','S-GLYCOSYLATION','polyphen_nobs','polyphen_normasa','polyphen_dvol','polyphen_dprop','polyphen_bfact','polyphen_hbonds','polyphen_avenhet','polyphen_mindhet','polyphen_avenint','polyphen_mindint','polyphen_avensit','polyphen_mindsit','polyphen_idpmax','polyphen_idpsnp','polyphen_idqmin','MOD_RES','REGION','INTERACTION_REGION','REQUIRED_FOR_INTER','Dst2Splice','Grantham','SpliceAI-acc-gain','SpliceAI-acc-loss','SpliceAI-don-gain','SpliceAI-don-loss','MMSp_acceptorIntron','MMSp_acceptor','MMSp_exon','MMSp_donor','MMSp_donorIntron','dbscSNV-ada_score','dbscSNV-rf_score','Charge','Volume','Hydrophobicity','Polarity','Ex','PAM250','BLOSUM','JM','HGMD2003','VB','Transition','COSMIC','COSMICvsSWISSPROT','HAPMAP','COSMICvsHAPMAP','ATP_binding_gbind','Ca2+_binding_gbind','DNA_binding_gbind','HEME_binding_gbind','Mg2+_binding_gbind','Mn2+_binding_gbind','RNA_binding_gbind']
MEDIAN_FEATURES = ['Conservation','TSSDistance','Eigen_PC_raw_coding','Eigen_raw_coding','GERPplus_plus_NR','GERPplus_plus_RS','GM12878_confidence_value','GM12878_fitCons_score','GenoCanyon_score','H1_hESC_confidence_value','H1_hESC_fitCons_score','HUVEC_confidence_value','HUVEC_fitCons_score','LINSIGHT','LIST_S2_score','LRT_Omega','LRT_score','MPC_score','MutationAssessor_score','SiPhy_29way_logOdds','bStatistic_x','integrated_confidence_value','integrated_fitCons_score','phastCons100way_vertebrate','phastCons17way_primate','phastCons30way_mammalian','phyloP100way_vertebrate','phyloP17way_primate','phyloP30way_mammalian','GDI','MSC_95CI','Selective_pressure','Clarks_distance','CDS_len','Number_of_paralogs','denovo_Zscore','RVIS','Indispensability_score','NearestExonJB_distance','s_het','gtex_Adipose_-_Subcutaneous','gtex_Adipose_-_Visceral_(Omentum)','gtex_Adrenal_Gland','gtex_Artery_-_Aorta','gtex_Artery_-_Coronary','gtex_Artery_-_Tibial','gtex_Bladder','gtex_Brain_-_Amygdala','gtex_Brain_-_Anterior_cingulate_cortex_(BA24)','gtex_Brain_-_Caudate_(basal_ganglia)','gtex_Brain_-_Cerebellar_Hemisphere','gtex_Brain_-_Cerebellum','gtex_Brain_-_Cortex','gtex_Brain_-_Frontal_Cortex_(BA9)','gtex_Brain_-_Hippocampus','gtex_Brain_-_Hypothalamus','gtex_Brain_-_Nucleus_accumbens_(basal_ganglia)','gtex_Brain_-_Putamen_(basal_ganglia)','gtex_Brain_-_Spinal_cord_(cervical_c-1)','gtex_Brain_-_Substantia_nigra','gtex_Breast_-_Mammary_Tissue','gtex_Cells_-_Cultured_fibroblasts','gtex_Cells_-_EBV-transformed_lymphocytes','gtex_Cervix_-_Ectocervix','gtex_Cervix_-_Endocervix','gtex_Colon_-_Sigmoid','gtex_Colon_-_Transverse','gtex_Esophagus_-_Gastroesophageal_Junction','gtex_Esophagus_-_Mucosa','gtex_Esophagus_-_Muscularis','gtex_Fallopian_Tube','gtex_Heart_-_Atrial_Appendage','gtex_Heart_-_Left_Ventricle','gtex_Kidney_-_Cortex','gtex_Kidney_-_Medulla','gtex_Liver','gtex_Lung','gtex_Minor_Salivary_Gland','gtex_Muscle_-_Skeletal','gtex_Nerve_-_Tibial','gtex_Ovary','gtex_Pancreas','gtex_Pituitary','gtex_Prostate','gtex_Skin_-_Not_Sun_Exposed_(Suprapubic)','gtex_Skin_-_Sun_Exposed_(Lower_leg)','gtex_Small_Intestine_-_Terminal_Ileum','gtex_Spleen','gtex_Stomach','gtex_Testis','gtex_Thyroid','gtex_Uterus','gtex_Vagina','gtex_Whole_Blood','haplo','haplo_imputed','GC','CpG','motifECount','motifEHIPos','motifEScoreChng','minDistTSS','minDistTSE','priPhCons','mamPhCons','verPhCons','priPhyloP','mamPhyloP','verPhyloP','targetScan','mirSVR-Score','mirSVR-E','mirSVR-Aln','cHmm_E1','cHmm_E2','cHmm_E3','cHmm_E4','cHmm_E5','cHmm_E6','cHmm_E7','cHmm_E8','cHmm_E9','cHmm_E10','cHmm_E11','cHmm_E12','cHmm_E13','cHmm_E14','cHmm_E15','cHmm_E16','cHmm_E17','cHmm_E18','cHmm_E19','cHmm_E20','cHmm_E21','cHmm_E22','cHmm_E23','cHmm_E24','cHmm_E25','GerpRS','GerpRSpval','GerpN','GerpS','tOverlapMotifs','motifDist','EncodeH3K4me1-sum','EncodeH3K4me1-max','EncodeH3K4me2-sum','EncodeH3K4me2-max','EncodeH3K4me3-sum','EncodeH3K4me3-max','EncodeH3K9ac-sum','EncodeH3K9ac-max','EncodeH3K9me3-sum','EncodeH3K9me3-max','EncodeH3K27ac-sum','EncodeH3K27ac-max','EncodeH3K27me3-sum','EncodeH3K27me3-max','EncodeH3K36me3-sum','EncodeH3K36me3-max','EncodeH3K79me2-sum','EncodeH3K79me2-max','EncodeH4K20me1-sum','EncodeH4K20me1-max','EncodeH2AFZ-sum','EncodeH2AFZ-max','EncodeDNase-sum','EncodeDNase-max','EncodetotalRNA-sum','EncodetotalRNA-max','Dist2Mutation','Freq100bp','Rare100bp','Sngl100bp','Freq1000bp','Rare1000bp','Sngl1000bp','Freq10000bp','Rare10000bp','Sngl10000bp','RemapOverlapTF','RemapOverlapCL']

def run_mi(l, ns):
    print('started')
    y = ns.lab[l].tolist()
    warnings.filterwarnings("ignore", category=UserWarning)
    mi = mutual_info_classif(ns.df, y, copy=True, discrete_features=True, random_state=1)
    return l, mi

def select_gene_features(data, labels):
    mgr = Manager()
    ns = mgr.Namespace()
    ns.df = data
    ns.lab = labels
    with Pool(labels.shape[1]) as p:
        mis = p.map(partial(run_mi, ns=ns), labels.columns.tolist())
    out = {mi[0]: mi[1] for mi in mis}
    out = pd.DataFrame.from_dict(out)
    out['max'] = out.apply(lambda r: max(r.tolist()), axis=1)
    p = np.percentile(out['max'].tolist(), [99])
    out['keep'] = out['max'].apply(lambda v: v >= p[0])
    return out

def run_mi_lp(l, ns):
    print('started')
    y = np.zeros(len(ns.lab))
    y[ns.lab == l] = 1
    warnings.filterwarnings("ignore", category=UserWarning)
    mi = mutual_info_classif(ns.df, y, copy=True, discrete_features=True, random_state=1)
    return l, mi

def select_gene_features_lp(data, labels, n_jobs=40):
    label_powerset = {}
    reverse_combinations = []
    last_id = 0
    labels_transformed = []
    for labels_applied in labels.values:
        label_string = "".join(map(lambda x: str(int(x)), labels_applied))

        if label_string not in label_powerset:
            label_powerset[label_string] = last_id
            reverse_combinations.append(labels_applied)
            last_id += 1

        labels_transformed.append(label_powerset[label_string])
    labels_transformed = np.array(labels_transformed)
    topls = list(set(sorted(labels_transformed,key=labels_transformed.tolist().count,reverse=True)))[0:80]

    mgr = Manager()
    ns = mgr.Namespace()
    ns.df = data
    ns.lab = labels_transformed
    with Pool(n_jobs) as p:
        mis = p.map(partial(run_mi_lp, ns=ns), topls)
    out = {mi[0]: mi[1] for mi in mis}
    out = pd.DataFrame.from_dict(out)
    out['max'] = out.apply(lambda r: max(r.tolist()), axis=1)
    p = np.percentile(out['max'].tolist(), [99])
    out['keep'] = out['max'].apply(lambda v: v >= p[0])
    return out

class DropCorrelatedTransformer(BaseEstimator, TransformerMixin):
    def __init__(self, correlation_cutoff=.95):
        self.correlation_cutoff = correlation_cutoff
        self.to_drop = None

    def fit(self, X, y=None):
        X = pd.DataFrame(X)
        cor_matrix = X.corr().abs()
        upper_tri = cor_matrix.where(np.triu(np.ones(cor_matrix.shape),k=1).astype(np.bool))
        self.to_drop = [column for column in upper_tri.columns if any(upper_tri[column] >= self.correlation_cutoff)]
        return self

    def transform(self, X, y=None):
        X = pd.DataFrame(X)
        X = X.drop(columns=self.to_drop)
        return X.to_numpy()

class SubsetTransformer(BaseEstimator, TransformerMixin):
    def __init__(self, max_columns):
        self.max_columns = max_columns
        self.keep_cols = None

    def fit(self, X, y=None):
        print(X.shape)
        cols = random.sample(list(range(X.shape[1])), round(X.shape[1] * self.max_columns))
        self.keep_cols = cols
        return self

    def transform(self, X, y=None):
        X = pd.DataFrame(X)
        X = X[self.keep_cols]
        print(X.shape)
        return X.to_numpy()

class FillNA(BaseEstimator, TransformerMixin):
    def __init__(self, fill=0):
        self.fill = 0

    def fit(self, X, y=None):
        self.X = X
        return self

    def transform(self, X, y=None):
        self.X = X
        return X.fillna(0)
    
    def get_feature_names_out(self, out):
        return self.X.columns.tolist()

class FastImputer(BaseEstimator, TransformerMixin):
    def __init__(self):
        pass

    def fit(self, X, y=None):
        print('fitting')
        self.medians = np.nanmean(X, axis=0)
        print('done fitting')
        return self

    def transform(self, X, y=None):
        X = X.copy()
        for i in tqdm(range(X.shape[1]), total = X.shape[1]):
            X[:,i][X[:,i] == np.nan] = self.medians[i]
        return X

def drop_allnan(data):
    for col in data.columns:
        if data[col].isna().sum() == len(data):
            data = data.drop(columns=col)
    return data

class BooltoInt(BaseEstimator, TransformerMixin):
    def __init__(self):
        pass

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        return X.replace([False, True], [0, 1]) 

def generate_preprocessor(numeric_features, categorical_features, bool_features, N_JOBS, do_categorical=True):

    categorical_transformer = OrdinalEncoder(handle_unknown='use_encoded_value', unknown_value=-1, encoded_missing_value=-1)

    numeric_transformer = Pipeline(steps=[
        ('imputer', FillNA()),
        ('scaler', MinMaxScaler(feature_range =(0, 1), clip=True))])

    transformers = [('numeric', numeric_transformer, numeric_features),
                    ('bool', BooltoInt(), bool_features)]
    if do_categorical:
        transformers.append(('cat', categorical_transformer, categorical_features))
    else:
        transformers.append(('remainder', 'passthrough'))

    preprocessor = ColumnTransformer(
        transformers=[
            ('numeric', numeric_transformer, numeric_features),
            ('cat', categorical_transformer, categorical_features),
        ])

    steps = [('initial', preprocessor)]
    if do_categorical:
        steps.append(('variance_threshold', VarianceThreshold(threshold=0)))
    preprocessor = Pipeline(steps=steps)
    return preprocessor

def preprocess(preprocessor, train_data, train_labels, test_data, quiet=False):
    for k, v in preprocessor.steps:
        if k == 'initial':
            start = time.time()
            v.fit(train_data)
            train_data = pd.DataFrame(v.transform(train_data), columns=v.get_feature_names_out())
            end = time.time()
            if not quiet:
                print(k + ' took ' + str(end - start) + ' to run.')
        elif k == 'oversampling' or k == 'undersampling':
            start = time.time()
            train_data, train_labels = v.fit_resample(train_data, train_labels)
            end = time.time()
            if not quiet:
                print(k + ' took ' + str(end - start) + ' to run.')
        else:
            start = time.time()
            train_data = v.fit_transform(train_data, train_labels)
            end = time.time()
            if not quiet:
                print(k + ' took ' + str(end - start) + ' to run')
    for k, v in preprocessor.steps:
        if k == 'initial':
            test_data = pd.DataFrame(v.transform(test_data), columns=v.get_feature_names_out())
        elif k == 'oversampling' or k == 'undersampling':
            continue
        else:
            test_data = v.transform(test_data)

    train_data = pd.DataFrame(train_data)
    test_data = pd.DataFrame(test_data)

    for col in train_data.columns:
        try:
            train_data[col] = train_data[col].astype('float')
            test_data[col] = test_data[col].astype('float')
        except:
            train_data[col] = train_data[col].astype('category')
            test_data[col] = test_data[col].astype('category')

    return train_data.to_numpy(), train_labels, test_data.to_numpy()


def transform(test_data, quiet=False):
    for k, v in preprocessor.steps:
        if k == 'initial':
            test_data = pd.DataFrame(v.transform(test_data), columns=v.get_feature_names_out())
        elif k == 'oversampling' or k == 'undersampling':
            continue
        else:
            test_data = v.transform(test_data)
    test_data = pd.DataFrame(test_data)

    for col in test_data.columns:
        try:
            test_data[col] = test_data[col].astype('float')
        except:
            test_data[col] = test_data[col].astype('category')

    return test_data.to_numpy()

def get_splits(seed, data, groups):
    random.Random(seed).shuffle(groups)
    fold_len = int(len(groups) * .2)
    fold_genes = np.array([groups[fold_len*i:fold_len*(i+1)] for i in range(5)])
    
    distance = 0
    for i in range(5):
        train_data = data[['cluster']].loc[data['cluster'].isin([g for fold in fold_genes[[f for f in range(5) if f != i]] for g in fold])]
        test_data = data[['cluster']].loc[data['cluster'].isin(fold_genes[i])]
        
        yield train_data.index.tolist(), test_data.index.tolist() 

def gbdtmo_trial(trial, n_classes=24, N_JOBS=40):
    GBDTMO_LIB = load_lib("/sc/arion/work/steind16/bin/GBDTMO/build/gbdtmo.so")
    out_dim = n_classes
    max_depth = trial.suggest_int('max_depth', 1, 50)
    max_leaves = trial.suggest_int('max_leaves', 2, 3002, step=20)
    min_samples = trial.suggest_int("min_samples", 5, 100)
    subsample = trial.suggest_float('subsample', .4, 1)
    lr = trial.suggest_float("lr", 0.01, 0.3)

    params = {
        'max_depth': max_depth,
        'max_leaves': max_leaves,
        'min_samples': min_samples,
        'subsample': subsample,
        'lr': lr,
        'num_threads': N_JOBS
    }
    print(params)

    return GBDTMulti(GBDTMO_LIB, out_dim=out_dim, params=params)

def boomer_trial(trial):

    max_rules = trial.suggest_int('max_rules', 200 , 2000, step=200)
    do_label_sampling = trial.suggest_categorical('do_label_sampling', [False, True])
    label_sampling = 'none'
    if do_label_sampling:
        label_sampling = 'without-replacement{num_samples=10}'
    feature_binning = trial.suggest_categorical('feature_binning', ['equal-width', 'equal-frequency'])
    shrinkage = trial.suggest_float('shrinkage', .1, 1, step=.1)
    loss = trial.suggest_categorical('loss', ['logistic-label-wise','logistic-example-wise','squared-error-label-wise','squared-hinge-label-wise'])

    return Boomer(max_rules=max_rules, label_sampling=label_sampling, feature_binning=feature_binning, shrinkage=shrinkage, loss=loss)

def xgboost_trial(trial, N_JOBS=40):
    # n_estimators = trial.suggest_int('n_estimators', 100, 1000, step=50)
    # learning_rate = trial.suggest_float("learning_rate", .001, .401, step=.001)
    min_child_weight = trial.suggest_int("min_child_weight", 0, 200)
    max_depth = trial.suggest_int('max_depth', 1, 200)
    colsample_bytree = trial.suggest_float('colsample_bytree', .4, 1, step=.1)
    colsample_bylevel = trial.suggest_float('colsample_bylevel', .4, 1, step=.1)
    max_delta_step = trial.suggest_int("max_delta_step", 0, 10)
    gamma = trial.suggest_float("gamma", 0, 100)
    reg_lambda = trial.suggest_int("reg_lambda", 0, 100)
    reg_alpha = trial.suggest_int("reg_alpha", 0, 100)
    subsample = trial.suggest_float("subsample", .4, 1, step=.1)

    params = {
        'n_estimators': 10,#n_estimators,
        # "learning_rate": learning_rate,
        "min_child_weight": min_child_weight,
        'max_leaves': 100,
        'max_depth': 5,# max_depth,
        'colsample_bytree': colsample_bytree,
        'colsample_bylevel': colsample_bylevel,
        "max_delta_step": max_delta_step,
        "gamma": gamma,
        "reg_lambda": reg_lambda,
        "reg_alpha": reg_alpha,
        "subsample": subsample,
        # "tree_method": 'gpu_hist',
        # "gpu_id": 0
    }

    model = XGBClassifier(n_jobs=N_JOBS, objective='multi:softmax', use_label_encoder=False)
    model.set_params(**params)
    return model

def lightgbm_trial(trial, objective, n_classes=None, N_JOBS=40):
    n_estimators = trial.suggest_int('n_estimators', 100, 1000, step = 50)
    num_leaves = trial.suggest_int('num_leaves', 2, 3002, step=20)
    learning_rate = trial.suggest_float("learning_rate", 0.01, 0.3)
    min_child_weight = trial.suggest_float("min_child_weight", 0.01, 20)
    min_child_samples = trial.suggest_int("min_child_samples", 5, 100)
    max_depth = trial.suggest_int('max_depth', -1, 20)
    colsample_bytree = trial.suggest_float('colsample_bytree', .4, 1)
    subsample = trial.suggest_float('subsample', .4, 1)
    subsample_freq = 0
    if subsample != 1:
        subsample_freq = 1
    reg_lambda = trial.suggest_int("reg_lambda", 0, 100)
    reg_alpha = trial.suggest_int("reg_alpha", 0, 100)
    is_unbalance = trial.suggest_int('is_unbalance', 0, 1)

    params = {
        'n_estimators': n_estimators,
        'num_leaves': num_leaves,
        'learning_rate': learning_rate,
        'min_child_weight': min_child_weight,
        'min_child_samples': min_child_samples,
        'max_depth': max_depth,
        'colsample_bytree': colsample_bytree,
        'subsample': subsample,
        'subsample_freq': subsample_freq,
        'reg_lambda': reg_lambda,
        'reg_alpha': reg_alpha,
        'is_unbalance': True if is_unbalance else False,
        'n_jobs': N_JOBS,
        'objective': objective,
        'verbose': -1   
    }

    if n_classes:
        params["num_class"] = n_classes

    classifier = LGBMClassifier()
    classifier.set_params(**params)
    return classifier

def sketchboost_trial(trial):
    ntrees = trial.suggest_int('ntrees', 100, 1000, step = 50)
    colsample = trial.suggest_float('colsample', .4, 1, step = .01)
    subsample = trial.suggest_float('subsample', .4, 1, step = .01)
    # max_depth = trial.suggest_int('max_depth', 2, 20)
    # lr = trial.suggest_float("lr", 0.01, 0.3, step = .01)
    lambda_l2 = trial.suggest_int("lambda_l2", 0, 100)
    use_hess = [False, True][trial.suggest_int("use_hess", 0, 1)]
    sketch_method = trial.suggest_categorical('sketch_method', ['filter', 'topk', 'rand', 'proj', None])

    return SketchBoost(
        'bce',
        ntrees=ntrees, 
        # lr=lr, 
        verbose=100, 
        lambda_l2=lambda_l2, 
        subsample=subsample, 
        colsample=colsample,
        # max_depth=max_depth, 
        use_hess=use_hess, 
        sketch_method=sketch_method
    )

def naivebayes_trial(trial):

    classifier = GaussianNB()
    return classifier

def LLSF_weight_mat(X,Y,optmParameter):
    
    #Optimal Parameters
    
        alpha            = optmParameter['alpha']
        beta             = optmParameter['beta']
        gamma            = optmParameter['gamma']
        maxIter          = optmParameter['maxIter']
        miniLossMargin   = optmParameter['minimumLossMargin']
    
    #Initialisation
    
        num_dim          = X.shape[1]
        XTX              = np.matmul(X.T,X)
        XTY              = np.matmul(X.T,Y)
        W_s              = np.matmul(inv(XTX + gamma * np.eye(num_dim)),XTY)
        W_s_1            = W_s
        eps              = 10**-8
        R                = cossim(np.transpose(Y + eps),np.transpose(Y + eps))
        Lip              = np.sqrt(2*np.power((LA.norm(XTX,'fro')),2) + np.power((LA.norm(alpha * R,'fro')),2))
        bk               = 1
        bk_1             = 1  
    
    #The Soft-thresholding function
    
        def softthres(W_t,lambd):
            W = np.maximum((W_t-lambd),lambd) - np.maximum(-W_t-lambd,lambd)
            return W
        
    #LLSF algorithm using Accelersted Proximal Gradient
        oldloss             = 0
        iteration           = 0
        while iteration <= maxIter:
            W_s_k           = W_s + ((bk_1 - 1)/bk) * (W_s - W_s_1)
            Gw_s_k          = W_s_k - ((1/Lip) * ((np.matmul(XTX,W_s_k) - XTY + alpha * np.matmul(W_s_k,R))))
            bk_1            = bk
            bk              = (1 + np.sqrt(4*bk**2 + 1))/2
            W_s_1           = W_s
            W_s             = softthres(Gw_s_k,beta/Lip)
            predictionLoss  = LA.norm((X@W_s - Y),'fro')
            correlation     = np.trace(np.matmul(R,np.matmul(W_s.T,W_s)))
            sparsity        = 1.0 - ( count_nonzero(W_s) / float(W_s.size) )  #sum(sum(W_s!=0))
            totalloss       = predictionLoss + alpha*correlation + beta*sparsity
            if np.absolute(oldloss - totalloss) <= miniLossMargin:
                break
            elif totalloss <=0:
                break
            else:
                oldloss = totalloss
            iteration+=1
        return W_s

class TimeoutException(Exception):   # Custom exception class
    pass

def break_after(seconds=2):
    def timeout_handler(signum, frame):   # Custom signal handler
        raise TimeoutException
    def function(function):
        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(seconds)
            try:
                res = function(*args, **kwargs)
                signal.alarm(0)      # Clear alarm
                return res
            except TimeoutException:
                print(u'Oops, timeout: %s sec reached.' % seconds, function.__name__)
            return
        return wrapper
    return function

def memory_limit(max_memory_gb, check_interval=0.1):
    max_memory_bytes = max_memory_gb * (1024 ** 3)

    def decorator(func):
        def wrapper(*args, **kwargs):
            p = psutil.Process()

            def check_memory():
                while p.is_running():
                    memory_usage = p.memory_info().rss
                    if memory_usage > max_memory_bytes:
                        raise MemoryError("Memory limit exceeded")
                    time.sleep(check_interval)

            # Start the memory check process
            memory_check_process = multiprocessing.Process(target=check_memory)
            memory_check_process.start()

            # Call the function
            result = func(*args, **kwargs)

            # Stop the memory check process after the function completes
            memory_check_process.terminate()
            memory_check_process.join()

            return result

        return wrapper

    return decorator

def calc_stats(y):
    num_instances = y.shape[0]
    num_labels = y.shape[1]

    c1, c0 = calcC0C1(y)
    max_C1 = max(c1)
    IRLbls = np.zeros(num_labels)
    ImRs = np.zeros(num_labels)
    SCUMBLEs = np.zeros(num_instances)

    mean_IR = 0;
    mean_ImR = 0;
    CVIR = 0;
    CVImR = 0;
    SCUMBLE = 0;

    mc = 0
    for j in range(num_labels):
        if c1[j] > mc:
            mc = c1[j]

    for j in range(num_labels):
        min_class_count = 0
        max_class_count = 0
        min_class_count = c1[j] if c1[j] < c0[j] else c0[j] 
        max_class_count = c1[j] if c1[j] > c0[j] else c0[j] 
        IRLbls[j] = mc / c1[j]
        ImRs[j] = max_class_count / min_class_count

    max_IRLb = max([v for i, v in enumerate(IRLbls) if c1[i] != 0])
    max_ImR = max([v for i, v in enumerate(ImRs) if min(c1[i], c0[i]) != 0])
    min_IRLb = min([v for i, v in enumerate(IRLbls) if c1[i] != 0])
    min_ImR = min([v for i, v in enumerate(ImRs) if min(c1[i], c0[i]) != 0])	
    cIR = sum([c1[i] != 0 for i in range(num_labels)])
    cImR = sum([min(c1[i], c0[i]) != 0 for i in range(num_labels)])
    mean_IR = sum([v for i, v in enumerate(IRLbls) if c1[i] != 0]) / cIR
    mean_ImR = sum([v for i, v in enumerate(ImRs) if min(c1[i], c0[i]) != 0]) / cImR

    CVIR = sum([(v - mean_IR) ** 2 for v in IRLbls if v != 0 and v < np.Inf])
    CVIR = math.sqrt(CVIR / (cIR - 1)) / mean_IR
    CVImR = sum([(v - mean_ImR) ** 2 for v in ImRs if v != 0 and v < np.Inf])
    CVIR = math.sqrt(CVImR / (cImR - 1)) / mean_ImR

    for i in range(num_instances):
        instance = y[i, :]
        pro = 1
        ave = 0
        count = 0
        for j in range(num_labels):
            if instance[j] == 1:
                pro = pro * IRLbls[j]
                ave = ave + IRLbls[j]
                count = count + 1
        if count == 0:
            SCUMBLEs[i] = 0
        else:
            ave = ave / count
            SCUMBLEs[i] = 1 - (1 / ave) * (pro ** (1 / count))

        SCUMBLE = SCUMBLE + SCUMBLEs[i]
    SCUMBLE = SCUMBLE / num_instances

    npl = len([v for v in c1 if v == 0]) / num_labels
    return {'SCUMBLE': SCUMBLE, 'SCUMBLEs': SCUMBLEs, 'CVIR': CVIR,
            'CVImR': CVImR, "mean_IR": mean_IR, "mean_ImR": mean_ImR,
            'IRLbls': IRLbls, 'ImRs': ImRs, 'npl': npl}

def calcC0C1(y):
    num_instances = y.shape[0]
    num_labels = y.shape[1]

    c1 = np.zeros(num_labels)
    c0 = np.zeros(num_labels)
    for i in range(num_instances):
        instance = y[i, :]
        for j in range(num_labels):
            if instance[j] == 1:
                c1[j] = c1[j] + 1
            else:
                c0[j] = c0[j] + 1
    return c1, c0

def conditional_random_int(min_val, max_val):
    if (min_val > max_val):
        return -1
    if min_val == max_val:
        return min_val
    return random.randint(0, 2147483647) % (max_val - min_val + 1) + min_val


def ml_random_oversample(y, p):

    sampled_indices = []
    num_instances = y.shape[0]
    num_labels = y.shape[1]
    num_add_instances = int(p * num_instances);
    imbalance_stats = calc_stats(y)
    c1, c0 = calcC0C1(y)

    low_labels = [imbalance_stats['IRLbls'][i] > imbalance_stats['mean_IR'] for i in range(num_labels)]
    min_bags = [[] for _ in range(num_labels)]

    for i in range(num_instances):
        instance = y[i, :]
        for j in range(num_labels):
            if instance[j] == 1 and low_labels[j]:
                min_bags[j].append(i)

    mc = -np.Inf;
    for c in c1:
        if mc < c:
            mc = c
		
    all_majority = True
    while num_add_instances > 0:
        
        all_majority = False if sum(low_labels) > 0 else True
        if all_majority:
            break
        
        for i in range(num_labels):
            if not low_labels[i]:
                continue
            
            new_IRLbl = mc / c1[i]
            if new_IRLbl <= imbalance_stats['mean_IR']:
                low_labels[i] = False
                continue

            k = conditional_random_int(0, len(min_bags[i]) - 1)
            sampled_indices.append(k)

            instance = y[k, :]
            for j in range(num_labels):
                if instance[j] == 1:
                    c1[j] = c1[j] + 1

            num_add_instances = num_add_instances - 1

    return sampled_indices

def remedial(y):
    y = y.copy()
    sampled_indices = []
    num_instances = y.shape[0]
    num_labels = y.shape[1]
    imbalance_stats = calc_stats(y)

    for i in range(num_instances):
        if imbalance_stats['SCUMBLEs'][i] > imbalance_stats['SCUMBLE']:
            for j in range(num_labels):
                if imbalance_stats['IRLbls'][j] <= imbalance_stats['mean_IR']:
                    y[i, j] = 0
                else:
                    y[i, j] = 1
    
    return y
