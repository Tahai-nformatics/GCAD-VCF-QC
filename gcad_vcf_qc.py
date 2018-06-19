#!/usr/bin/env python3
#-*- coding: utf-8 -*-

from argparse import ArgumentParser
from pysam import VariantFile

#import configparser
import csv
from collections import namedtuple
import time
import cProfile

from utils.VariantAnnotation.vflags import calcVFlags
import config as cfg


def extractSubSets(fam):
    """
    @return Samples Dict per Subset
    """
    SampleFamDetail = namedtuple('SampleFamDetail',['SampID','FID','SubjID','FA','MO','SEX','AFF','AD','AGE','ADSPWGS','Subset','Subgroup','Race_Ethnicity','SeqCtr','ExcludeFromZHet'])

    samples = dict()
    for sm in map(SampleFamDetail._make, csv.reader(open(fam, 'r'),delimiter='\t')):
        if sm.Subset in samples:
            samples[sm.Subset].add(sm.SampID)
        else:
            samples[sm.Subset]=set()
            samples[sm.Subset].add(sm.SampID)
    return samples



def main():
    argparser = ArgumentParser()
    grp_file_paths = argparser.add_argument_group(title='File paths')
    grp_file_paths.add_argument('--vcf', type=str, help='input VCF file', required=True)
    grp_file_paths.add_argument('--fam', type=str, help='fam file (with headers)', required=True)
    #grp_file_paths.add_argument('--outfile', type=str, help='filepath for output file', required=True)

    grp_settings = argparser.add_argument_group(title='Theshold Settings')
    grp_settings.add_argument('--min_dp', type=int, help='minimum depth (DP)', default=10, required=False)
    grp_settings.add_argument('--min_gq', type=int, help='maximum genotype quality (GQ)', default=20, required=False)
    grp_settings.add_argument('--isFam', type=int, help='is this family data, 1=yes, 0=no', default=1, required=False)
    grp_settings.add_argument('--minTranche', type=int, help='Minimum Tranche Score (Filter from VQSR)', default=99.7, required=False)
    grp_settings.add_argument('--miss_rate', type=int, help='max missingness threshold', default=0.2, required=False)
    grp_settings.add_argument('--max_dp', type=int, help='max DP threshold', default=500, required=False)
    grp_settings.add_argument('--hetz_lim1', type=int, help='', default=5, required=False)
    grp_settings.add_argument('--hetz_lim2', type=int, help='', default=6.101825, required=False)
    grp_settings.add_argument('--hwe_pval', type=int, help='', default=5e-06, required=False)
    grp_settings.add_argument('--hwe_maf', type=int, help='MAF threshold', default=0.01, required=False)

    args, extr = argparser.parse_known_args()

    cfg.MINDP = args.min_dp
    cfg.MINGQ = args.min_gq
    cfg.minTranche = args.minTranche
    cfg.miss_rate = args.miss_rate
    cfg.max_dp = args.max_dp
    cfg.isFam = args.isFam
    cfg.hetz_lim1 = args.hetz_lim1
    cfg.hetz_lim2 = args.hetz_lim2
    cfg.hwe_pval = args.hwe_pval
    cfg.hwe_maf = args.hwe_maf

    samplesDict = extractSubSets(args.fam) # returns dict
    #{k:rec.samples[k] for k in samples['sub1']}
    #{k:rec.samples[k] for k in samples['sub1'] if k in rec.samples}
    #{key: d[key] for key in d.viewkeys() & l}

    vcf_in = VariantFile(args.vcf)
    #vcf_out = VariantFile(args.outfile, 'w', header=vcf_in.header, threads=4)


    # Save groups of samples as an intersecting set
    for k,v in samplesDict.items():
        #extraction_set.append( list(set(vcf_in.header.samples) & set(v) ))
        #samplesDict[k] = list(set(vcf_in.header.samples) & set(v) )

        #method 3
        #samplesDict[k] = {'list': list(set(vcf_in.header.samples) & set(v) ),'dict':dict() }

        #method 4
        #samplesDict[k] = set(vcf_in.header.samples) & set(v)

        #method 5
        samplesDict[k] = {'set': set(vcf_in.header.samples) & v,'dict':dict() }

    ct=0
    start = time.time()
    for rec in vcf_in.fetch():
        #vcf_out.write(rec)
        if ct>99:break

        print(str(rec.pos) + '\t', end='')

        # Method 1
        #for sm_list in extraction_set:
        #    vf = calcVFlags1(rec.samples, rec.filter, sm_list)
        #    print("VFLAGS_{}={};".format('', vf), end='')

        # Method 2
        #for subset, sm_list in samplesDict.items():

            #gss = dict() # slice rec samples
            #for key,sm in rec.samples.items():
                #if key not in sm_list: continue ## !SLOW!
                #gss[key] = sm

            #vf = calcVFlags(gss, rec.filter)

         #   print("VFLAGS_{}={};".format(subset,vf), end='')

        #print()

        # Method 3
        #for key,sm in rec.samples.items():
            #for subset, sm_list in samplesDict.items():
                #if key in sm_list['list']: ## !SLOW!
                    #samplesDict[subset]['dict'][key] = sm

        #for subset, sm_list in samplesDict.items():
            #vf = calcVFlags(sm_list['dict'], rec.filter)
            #print("VFLAGS_{}={};".format(subset, vf), end='')

        # Method 4 - set()
        #for subset, sm_set in samplesDict.items():

            #gss = dict() # slice rec samples
            #for key,sm in rec.samples.items():

                #if key in sm_set:
                    #gss[key] = sm

            #vf = calcVFlags(gss, rec.filter)
            #print("VFLAGS_{}={};".format(subset, vf), end='')

        # Method 5 - set, one pass
        for key,sm in rec.samples.items():
            for subset, sm_set in samplesDict.items():
                if key in sm_set['set']:
                    samplesDict[subset]['dict'][key] = sm

        for subset, sm_list in samplesDict.items():
            vf = calcVFlags(sm_list['dict'], rec.filter)
            print("VFLAGS_{}={};".format(subset, vf), end='')


        print()
        ct += 1

    end = time.time()
    print(end - start)


if __name__ == "__main__":
    main()
    #cProfile.run('main()', None, 'cumtime')
