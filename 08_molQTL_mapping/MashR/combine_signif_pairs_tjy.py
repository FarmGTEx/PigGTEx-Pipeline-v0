#!/usr/bin/env python3
from multiprocessing import Pool
import pandas as pd
import numpy as np
import argparse
import os
import gzip


def read_dataframe(file):
    df = pd.read_csv(file, sep='\t', usecols=['variant_id', 'phenotype_id'],memory_map=True)
    return df
 
def read_file(split_file_path):
    file_list = split_file_path
    thread_num = len(file_list)
    if thread_num > 24:
        thread_num = 24
    p = Pool(thread_num)
    res = p.map(read_dataframe, file_list)
    p.close()
    p.join()
    df = pd.concat(res, axis=0)
    return df


parser = argparse.ArgumentParser(description='Combine significant pairs from multiple eQTL mapping runs.')
parser.add_argument('signifpair_list_file', help="File listing of 'signifpairs' outputs from eQTL pipeline.")
parser.add_argument('prefix', help='Prefix for output file: <prefix>.combined_signifpairs.txt.gz')
# parser.add_argument('--chunks', default=None, type=int, help='')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

with open(args.signifpair_list_file) as f:
    file_paths = f.read().strip().split('\n')

print('Loading significant pairs.')
chunks = int(np.ceil(len(file_paths) / 24))
dfs = []
if chunks is None:
    for f in file_paths:
        print(f)
        dfs.append(pd.read_csv(f, sep='\t', usecols=['variant_id', 'phenotype_id'],memory_map=True))
    dfs = pd.concat(dfs, axis=0)
    dfs.drop_duplicates(inplace=True)
else:
    chunk_size = int(np.ceil(len(file_paths) / chunks))
    for i in np.arange(chunks):
        print('  * loading chunk {}/{}'.format(i+1, chunks), flush=True)
        fp = file_paths[i*chunk_size:(i+1)*chunk_size]
        print(fp)
        dfchunk = read_file(fp)
        dfchunk.drop_duplicates(inplace=True)
        dfs.append(dfchunk)
        if len(dfs) > 1:
            tmp = pd.concat(dfs, axis=0)
            tmp.drop_duplicates(inplace=True)
            dfs = []
            dfs.append(tmp)
    dfs = dfs[0]

print('Sorting significant pairs.')
dfs['chr'] = dfs['variant_id'].apply(lambda x: x.split('_',1)[0])
dfs['pos'] = dfs['variant_id'].apply(lambda x: int(x.split('_',2)[1]))
dfs = dfs.sort_values(['chr', 'pos', 'phenotype_id'])

print('Writing output.')
with gzip.open(os.path.join(args.output_dir, args.prefix+'.combined_signifpairs.txt.gz'), 'wt', compresslevel=6) as f:
    dfs.to_csv(f, sep='\t', index=False)
    
    