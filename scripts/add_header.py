import pandas as pd
import sys

datafile = sys.argv[1]
headerfile = sys.argv[2]
outfile = sys.argv[3]

csqheader = pd.read_csv(headerfile, header=None, names=["_","name"], sep="\t")["name"].to_list()

csqheaderdedup = []
for i in range(len(csqheader)):
    if csqheader[i] in ['APPRIS', 'TSL']:
        csqheaderdedup.append(csqheader[i]+str(i))
    else:
        csqheaderdedup.append(csqheader[i])     
  
vcfheader = ["CHROM","POS","ID","REF","ALT"]
vcfheader+=csqheaderdedup

data = pd.read_csv(datafile, sep='\t', header=None)
data.columns = vcfheader
data.to_csv(outfile, index=None)
