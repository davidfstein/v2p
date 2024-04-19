import pandas as pd
import sys
from tqdm import tqdm
tqdm.pandas()

infile = sys.argv[1]
outfile = sys.argv[2]
gene = sys.argv[3] if len(sys.argv) >= 4 else ''

data = pd.read_csv(infile, low_memory=False)
if gene:
    data = data.loc[data['SYMBOL'] == gene].reset_index(drop=True)
data = data.drop_duplicates(subset=['ID'])

##take the ones with protein info
mydata=data[(data['BIOTYPE']=="protein_coding") & (data['Protein_position'].notnull()) & (data['Protein_position']!=".")].copy()

##Relative position: pos/len*10 (gives number btw 0-10)
def calc_pos(x):
    if x=="-":
        return None
    else:
        try:
            splitted=x.split("/")
            pos_info=splitted[0]
            prot_len=splitted[1]
        except:
            raise ValueError(x, splitted)
        positions=pos_info.split("-")
        if len(positions)==1: #single position, not a range
            if int(pos_info) > int(prot_len): ## if stop-gain: 503/502
                return 1

            else:
                return float(pos_info)/float(prot_len) 

        elif len(positions)==2: ## if position is a range
            start_pos=positions[0]
            end_pos=positions[1]
            if start_pos=="?" or start_pos=="":
                if end_pos!='':
                    return float(end_pos)/float(prot_len)
            elif int(start_pos) > int(prot_len):
                return 1
            else:
                return float(start_pos)/float(prot_len)

### Calculate relative positions (0-10): pos/len*10
mydata["rel_cDNA_pos"]=mydata["cDNA_position"].progress_apply(calc_pos)
mydata["rel_CDS_pos"]=mydata["CDS_position"].progress_apply(calc_pos)
mydata["rel_prot_pos"]=mydata["Protein_position"].progress_apply(calc_pos)

##Relative position: pos/len*10 (gives number btw 0-10)
def aa_pos(x):
    if x=="-":
        return "NULL"
    else:
        splitted=x.split("/")
        pos_info=splitted[0]
        prot_len=splitted[1]
        positions=pos_info.split("-")
        if len(positions)==1: #single position, not a range
            if int(pos_info) > int(prot_len): ## if stop-gain: 503/502
                return int(prot_len)
            else:
                return int(pos_info)

        elif len(positions)==2: ## if position is a range
            start_pos=positions[0]
            end_pos=positions[1]
            if start_pos=="?" or start_pos=="":
                if end_pos!='':
                    return int(end_pos)
            elif int(start_pos) > int(prot_len):
                return int(prot_len)
            else:
                return int(start_pos)

mydata["AA_pos"]=mydata["Protein_position"].progress_apply(aa_pos)


# ### Ref and Alt aa types
def ref_aa(x):
    if x=="-" or pd.isna(x):
        return "NULL"
    else:
        splitted=x.split("/")
        if len(splitted)==1: #synonymous
            return str(splitted[0])
        elif len(splitted)==2:
            return str(splitted[0])[0]
mydata["Ref_aa"]=mydata["Amino_acids"].progress_apply(ref_aa)

def alt_aa(x):
    if x=="-" or pd.isna(x):
        return "NULL"
    else:
        splitted=x.split("/")
        if len(splitted)==1: #synonymous
            return str(splitted[0])
        elif len(splitted)==2:
            if "*" in splitted[1]:
                return "*"
            elif "X" in splitted[1]:
                return "X"
            else:
                return str(splitted[1])[0]
mydata["Alt_aa"]=mydata["Amino_acids"].progress_apply(alt_aa)

mydata = mydata[["ID","rel_cDNA_pos","rel_CDS_pos","rel_prot_pos","AA_pos","Ref_aa","Alt_aa"]]
out = data.merge(mydata, on='ID', how='left')
out = out.replace(str('.'), None)
out.to_parquet(outfile, index=False)

