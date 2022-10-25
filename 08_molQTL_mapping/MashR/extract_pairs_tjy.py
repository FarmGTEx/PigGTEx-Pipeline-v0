#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import os
import gzip
import feather

parser = argparse.ArgumentParser(description='Extract variant-gene pairs from list of associations.')
parser.add_argument('allpairs_list_file', help="File listing of 'signifpairs' outputs from eQTL pipeline.")
# parser.add_argument('input_pairs', help="'Nominal' output from FastQTL.")
parser.add_argument('extract_pairs', help="File containing list variant-gene pairs to extract.")
parser.add_argument('prefix', help='Prefix for output file: <prefix>.extracted_pairs.{txt.gz|ft}')
parser.add_argument('--chrom', default=None, type=int, help='Chromosome')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
parser.add_argument('--feather', action='store_true', help='Write output in feather format')
args = parser.parse_args()


print('Loading list of pairs to extract.')
with open(args.allpairs_list_file) as f:
    file_paths = f.read().strip().split('\n')

print('Loading significant pairs.')
ref_pairs_df = pd.read_csv(args.extract_pairs, sep='\t', usecols=['variant_id', 'phenotype_id', 'chr'], dtype={'variant_id':str, 'phenotype_id':str, 'chr':int}, memory_map=True)
if args.chrom is not None:
    print(args.chrom)
    ref_pairs_df = ref_pairs_df[ref_pairs_df['chr']==args.chrom]
ref_pair_ids = set(ref_pairs_df['variant_id']+','+ref_pairs_df['phenotype_id'])

print('Parsing variant-gene pairs.')
signif_df = []
for f in file_paths:
    print(f)
    for i,chunk in enumerate(pd.read_csv(f, sep='\t', iterator=True, chunksize=1000000, index_col=1,
        dtype={'phenotype_id':str, 'variant_id':str, 'tss_distance':np.int32,
            'ma_samples':np.int32, 'ma_count':np.int32, 'af':np.float32,
            'pval_nominal':np.float64, 'slope':np.float32, 'slope_se':np.float32},memory_map=True)):

        ix = (chunk.index.values+','+chunk['phenotype_id']).apply(lambda x: x in ref_pair_ids)
        signif_df.append(chunk.iloc[ix.values])
        print('  * chunks processed: {0:d}'.format(i+1), end='\r', flush=True)

signif_df = pd.concat(signif_df, axis=0)
signif_df = signif_df.reset_index()
signif_df.index = signif_df['variant_id']+','+signif_df['phenotype_id']

print('\nSorting extracted pairs.')
output_df = pd.DataFrame(index=ref_pairs_df['variant_id']+','+ref_pairs_df['phenotype_id'])
output_df = output_df.join(signif_df[['pval_nominal', 'slope', 'slope_se']])
output_df.index.name = 'pair_id'

print('Writing output.')
if args.feather:
    output_df = output_df.reset_index()
    feather.write_dataframe(output_df, os.path.join(args.output_dir, args.prefix+'.extracted_pairs.ft'))
else:
    with gzip.open(os.path.join(args.output_dir, args.prefix+'.extracted_pairs.txt.gz'), 'wt', compresslevel=6) as f:
        output_df.to_csv(f, sep='\t', na_rep='NA', float_format='%.6g')
        

