import pandas as pd
import sys

vep_data = pd.read_parquet(sys.argv[1])
additional_data = pd.read_parquet(sys.argv[2])
out = vep_data.merge(additional_data, on='ID')
out.to_parquet(sys.argv[3], index=None)
