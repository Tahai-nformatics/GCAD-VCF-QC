#!/usr/bin/env python3
import os
import glob
from pathlib import Path
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='merge subset snv from intervals')

parser.add_argument('--dir',      help='input directory, recursive search', required=True)
parser.add_argument('--subset',   help='subset', required=True)
parser.add_argument('--outfile',  help='output filename', required=True)
args = parser.parse_args()

# Define the sorter
sorter = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
          'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
          'chr20','chr21','chr22',
         ]

def get_merged_csv(flist, **kwargs):
  """
   Read-in and concat files in glob
  """
  return pd.concat([pd.read_csv(f, sep='\t', **kwargs) for f in flist], sort=False)


fstr = 'summary.snv.*' + args.subset + '.tsv'
#fmask = os.path.join(args.dir, fstr)
#gg = glob.glob(fmask, recursive=True)
gg = list(Path(args.dir).rglob(fstr))
fcount = len(gg)

if fcount > 1:
  print("Merging {} files from {}, matching {}, into {}".format(len(gg), args.dir, fstr, args.outfile))
  df = get_merged_csv(gg)

  # Re-format floats
  for c in ['CallRate', 'MAF','AltAF', 'MeanDepth','propMI']:
     df[c] = df[c].map(lambda n: "{0:.6f}".format(n))

  df['ABHet'] = df['ABHet'].map(lambda n: "{0:.4f}".format(n) if not pd.isna(n) else 'NA' )
  df['QUAL'] = df['QUAL'].map(lambda n: "{0:.2f}".format(n) if not pd.isna(n) else 'NA' )

  # alter CHR column to allow sorting as in sorter list
  df.CHR = df.CHR.astype("category")
  df.CHR.cat.set_categories(sorter, inplace=True)

  # sort, save
  df.sort_values(['CHR','POS']).to_csv(args.outfile, index=False, sep='\t')
