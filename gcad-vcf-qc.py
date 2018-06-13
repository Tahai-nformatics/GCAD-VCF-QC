#!/usr/bin/env python3
#-*- coding: utf-8 -*-

from argparse import ArgumentParser
from pysam import VariantFile
from subprocess import check_output
#import configparser
import math

MINDP=0
MINGQ=0

def calc_InbreedingCoeff(aa, ab, bb):
    """
    Inbreeding Coefficient - Although the name Inbreeding Coefficient suggests 
    it is a measure of inbreeding, Inbreeding Coefficient measures the excess heterozygosity 
    at a variant site. It can be used as a proxy for poor mapping (sites that have 
    high Inbreeding Coefficients are typically locations in the genome where the mapping 
    is bad and reads that are in the region mismatch the region because they belong elsewhere). 
    At least 10 samples are required (preferably many more) in order for this annotation to be 
    calculated properly.

    We calculate Inbreeding Coefficient as
            # observed heterozygotes
        1 − ————————————————————————
            # expected heterozygotes

    # expected heterozygotes = 2pq * N
    """
    N = aa+ab+bb
    p = ( (2*aa)+ab )/ (2*N)
    q = 1 - p
    return 1 - ( ab / (2*p*q*N) )

def calc_ExcessHet(aa, ab, bb):
    """
    ExcessHet - Phred-scaled p-value for exact test of excess heterozygosity
    This annotation estimates excess heterozygosity in a population of samples. 
    It is related to but distinct from InbreedingCoeff, which estimates evidence 
    for inbreeding in a population. ExcessHet scales more reliably to large cohort sizes.

    p=freq of allele a in population
    q=1-p
    """
    N = aa+ab+bb

    result = '.'
    maf = 0
    if N:
        maf = ((2 * bb) + ab)/(2*N)
        #p =   ((2 * aa) + ab)/(2*N)
        #q = 1 - p

    if maf > 0.5: maf = 1 - maf

    if (ab == N):
        result = '.'
    elif ( (maf >0) and (ab >0) ):
        hetExpct = 2 * maf * ( 1 - maf)
        hetObs   = ab / N
        t  = hetObs - hetExpct
        chisq = (N * t * t) / ( hetObs * (1 - hetObs) )
        result = math.sqrt(chisq)
        if t < 0:
            result = -result

    return [result, maf]

def calc_ExcessHetOLD(aa, ab, bb):
    zhet, hetz_maf = (check_output(
                                ["perl",
                                    "-e",
                                    "use lib '/home/ottov/gcad-qc-pipeline';use ExcessHeterozygosity;@a=ExcessHeterozygosity::excess_het(%d,%d,%d);print join('\t',@a);" %
                                    (aa, ab, bb)])
                            .decode('ascii').split('\t'))
    return [zhet, hetz_maf]

def calc_ABHet(ad, dp):
    """
    AlleleBalanceHet - Allele Balance for heterozygous calls (ref/(ref+alt))
    {# REF reads from heterozygous samples}/{# REF + ALT reads from heterozygous samples}

    """
    if dp:
        return ad/dp
    else:
        return 0

def count_gt(samples):
   """
   count_gt -
   """
   obs_hts = 0  # zHet01
   obs_hom1 = 0 # zHet00
   obs_hom2 = 0 # zHet11

   for sm in samples.keys():
     if (samples[sm]['GT'] == (None, None)
         or
         samples[sm]['DP'] == None):
         continue
     elif (samples[sm]['DP'] < MINDP
         or samples[sm]['GQ'] < MINGQ
        ):
        continue

     if (samples[sm]['GT'] == (0,1)
       or samples[sm]['GT'] == (1,0)):
          obs_hts += 1
     elif samples[sm]['GT'] == (0,0):
          obs_hom1 += 1
     elif samples[sm]['GT'] == (1,1):
          obs_hom2 += 1
     #else: not counted

   return [obs_hom1, obs_hts, obs_hom2]

def main():
    argparser = ArgumentParser()
    grp_file_paths = argparser.add_argument_group(title='File paths')
    grp_file_paths.add_argument('--vcf', type=str, help='input VCF file', required=True)
    #grp_file_paths.add_argument('--outfile', type=str, help='filepath for output file', required=True)

    grp_settings = argparser.add_argument_group(title='Settings')
    grp_settings.add_argument('--min_dp', type=int, help='minimum depth (DP)', default=10, required=False)
    grp_settings.add_argument('--min_gq', type=int, help='maximum genotype quality (GQ)', default=20, required=False)
    args, extr = argparser.parse_known_args()

    MINDP = args.min_dp
    MINGQ = args.min_gq

    vcf_in = VariantFile(args.vcf)
    #vcf_out = VariantFile(args.outfile, 'w', header=vcf_in.header)

    for rec in vcf_in.fetch():
        #vcf_out.write(rec)
        [obs_hom1, obs_hets, obs_hom2] = count_gt(rec.samples)

        [zhet1, hetz_maf] = calc_ExcessHetOLD(obs_hom1, obs_hets, obs_hom2 )
        [zhet2, maf] = calc_ExcessHet(obs_hom1, obs_hets, obs_hom2 )



        if zhet1 != '.':zhet1="{:.4f}".format(float(zhet1))
        if zhet2 != '.':zhet2="{:.4f}".format(float(zhet2))

        if zhet1 == zhet2: continue

        print("zhet1:={};\tzhet2:={};".format(zhet1, zhet2))

if __name__ == "__main__":
    main()
