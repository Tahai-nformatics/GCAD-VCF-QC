#!/usr/bin/env python3

import pandas as pd
import argparse
from collections import OrderedDict
import subprocess

# dict for chr start position along linear chr
# hg38
cLookup = OrderedDict()
cLookup['chr1'] = 248956422
cLookup['chr2'] = 242193529
cLookup['chr3'] = 198295559
cLookup['chr4'] = 190214555
cLookup['chr5'] = 181538259
cLookup['chr6'] = 170805979
cLookup['chr7'] = 159345973
cLookup['chr8'] = 145138636
cLookup['chr9'] = 138394717
cLookup['chr10'] = 133797422
cLookup['chr11'] = 135086622
cLookup['chr12'] = 133275309
cLookup['chr13'] = 114364328
cLookup['chr14'] = 107043718
cLookup['chr15'] = 101991189
cLookup['chr16'] = 90338345
cLookup['chr17'] = 83257441
cLookup['chr18'] = 80373285
cLookup['chr19'] = 58617616
cLookup['chr20'] = 64444167
cLookup['chr21'] = 46709983
cLookup['chr22'] = 50818468

parser = argparse.ArgumentParser(description='Process a file to partitions')

parser.add_argument('csv', help='companion file')
parser.add_argument('vcf', help='vcf file')
parser.add_argument('chr', help='chromosome filter')
parser.add_argument('cuts', help='number of partitions')
args = parser.parse_args()


def pos_is_empty(pos):
   """
   pos_is_empty - use bcftools to check if the given position
                  is free of indels which would be broken by
                  the cut point
   """
   output = subprocess.check_output('bcftools view -V snps -H -r {CHR}:{REGION} {VCF} | wc -l'.format(CHR=args.chr,REGION=pos, VCF=args.vcf),                   shell=True)
   return int(output) == 0

def main():
 """

 """
 # determine bin for approx equal quantity bins (quantiles)
 g = pd.qcut(dat1['POS'], int(args.cuts), retbins=True, labels=False)

 # Convert results to int
 f = list(map(int, g[1] ))
 f[0]=1

 rd=[]

 for idx in range(0,len(f) - 1):

  # end_pos is next interval
  if idx < len(f) - 2:
      end_pos = f[idx+1] - 1

      ct = 0
      while not pos_is_empty(end_pos):
        end_pos += 1
        ct += 1

      if idx < len(f):
        f[idx+1] += ct

  else:
      end_pos = cLookup[chr]

  g = dict(
         chr   = chr,
         begin = f[idx],
         end   = end_pos,
          )
  rd.append(g)

 for im in rd:
   print ("%s:%s-%s" % (im['chr'],im['begin'],im['end']))


if __name__ == '__main__':
 chr=args.chr

 dat1 = pd.read_csv(args.csv, usecols=['POS'], sep='\t')

 main()
