import sys
import os
import pandas as pd

parfile=sys.argv[1]
outfile=sys.argv[2]
df = pd.read_parquet(parfile)
df.to_csv(outfile,compression="gzip",sep="\t",index=False)
