import pandas as pd
import numpy as np
import sys

vep_data = pd.read_parquet(sys.argv[1])
additional_data = None
try:
    additional_data = pd.read_parquet(sys.argv[2])
except:
    additional_data = pd.read_csv(sys.argv[2])
out = vep_data.merge(additional_data, on='ID')
try:
    out.to_parquet(sys.argv[3], index=None)
except:
    out.to_csv(sys.argv[3], index=None)
