#!/usr/bin/env python3
#-*- coding: utf-8 -*-

from argparse import ArgumentParser
from pysam import VariantFile

#import configparser
import csv
import os.path
from collections import namedtuple, OrderedDict
import time
import cProfile

from utils.VariantAnnotation.vflags import calcVA
import config as cfg

import warnings
warnings.simplefilter('always')

def extractSubSets(fam):
    """
    extractSubSets - Creates a OrderedDict() where Subset groupings are the key, 
                     values are SampID within it
    @return Samples Dict per Subset
    """
    SampleFamDetail = namedtuple('SampleFamDetail',['SampID','FID','SubjID','FA','MO','SEX','AFF','AD','AGE','ADSPWGS','Subset','Subgroup','Race_Ethnicity','SeqCtr','ExcludeFromZHet'])

    samples = OrderedDict()
    ct = 0
    with open(fam, 'r') as fam_file:
        for sm in map(SampleFamDetail._make, csv.reader(fam_file, delimiter='\t')):
            if sm.Subset in samples:
                samples[sm.Subset].add(sm.SampID)
            else:
                samples[sm.Subset] = set()
                samples[sm.Subset].add(sm.SampID)

            ct += 1
    return samples, ct

def write_subset_stats(subset, rec, vf, abhet, passing, failing, missing, gt_failed, depth_sum, clean_obs):
    """
    """
    outfile = 'summary.snv.{}_{}.out'.format(rec.contig, subset)
    newfile = not os.path.exists(outfile)


    with open(outfile, 'a') as csvfile:
        fieldnames = ['CHR','POS',
                      'Pass00','Pass01','Pass11',
                      'Fail00','Fail01','Fail11',
                      'Missing','GT_Failed',
                      'Clean00','Clean01','Clean11','Mono','CallRate','CallBad','GATKPass','MAF','AltAF',
                      'MeanDepth','HiDepth','ABHet','Mend_Incon','Mend_pairs','propMI','MultiAllele','FilteredOut',
                      'VFLAGS','rsID','RefAllele','AltAllele','QUAL','FILTER',
                      ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames , delimiter='\t')

        if newfile:
            writer.writeheader()

        # field calculations
        sum_clean = sum(clean_obs)

        # callrate
        callrate = 1 - (missing + gt_failed) / (missing + gt_failed + sum_clean)

        # MAF / AltAF; altAF = maf
        maf = 0
        temp = 2 * sum_clean

        if temp > 0:
            maf = (clean_obs[1] + 2 * clean_obs[2]) / temp
            if maf > 0.5:
                maf = 1 - maf

        maf = "{0:.6f}".format(maf)

        # MeanDepth
        mean_depth = depth_sum / sum_clean if sum_clean > 0 else 0

        writer.writerow({'CHR': rec.contig, 'POS': rec.pos,
                         'Pass00':passing[0],'Pass01':passing[1],'Pass11':passing[2],
                         'Fail00':failing[0],'Fail01':failing[1],'Fail11':failing[2],
                         'Missing':missing, 'GT_Failed':gt_failed,
                         'Clean00':clean_obs[0],'Clean01':clean_obs[1],'Clean11':clean_obs[2],
                         'Mono':int(3 in vf),
                         'CallRate':"{0:.6f}".format(callrate),
                         'CallBad': int(callrate <= (1 - cfg.miss_rate)),
                         'GATKPass':int(1 not in vf),
                         'MAF':maf,'AltAF':maf,
                         'MeanDepth':"{0:.6f}".format(mean_depth),'HiDepth':int(mean_depth > cfg.max_dp),
                         'ABHet':abhet,
                         'Mend_Incon':0,'Mend_pairs':0,'propMI':0,
                         'MultiAllele':0,'FilteredOut':0,
                         'VFLAGS':",".join(map(str,vf)),
                         'rsID':rec.id,'RefAllele':rec.ref,'AltAllele':rec.alts[0],
                         'QUAL':"{0:.2f}".format(rec.qual),'FILTER':",".join(rec.filter.keys()),
                         })


    return

def main():
    argparser = ArgumentParser()
    grp_file_paths = argparser.add_argument_group(title='File paths')
    grp_file_paths.add_argument('--vcf', type=str, help='input VCF file', required=True)
    grp_file_paths.add_argument('--fam', type=str, help='fam file (with headers)', required=True)
    #grp_file_paths.add_argument('--outfile', type=str, help='filepath for output file', required=True)

    grp_overrides = argparser.add_argument_group(title='bcf overrides')
    grp_overrides.add_argument('--region', type=str, help='restrict to VCF region e.g. chr3:100-200', default=None, required=False)

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

    createSampleAnnotation(args.fam)
    samplesDict, famCt = extractSubSets(args.fam) # returns dict
    print("[FAM] Found {} subsets: {}; for {} sampIDs".format(len(samplesDict.keys()), list(samplesDict.keys()), famCt))
    print("[FAM] {}".format([  "{}:{}".format(k, len(samplesDict[k]))  for k in samplesDict.keys()]))
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

    print("[VCF] contains {} samples".format(len(vcf_in.header.samples)))

    ct=0
    start = time.time()
    rChr = None
    rStart = None
    rEnd = None

    if args.region:
        rChr = args.region.split(':')[0]
        rStart = int(args.region.split(':')[1].split('-')[0]) - 1
        rEnd = int(args.region.split(':')[1].split('-')[1])

    for rec in vcf_in.fetch(rChr, rStart, rEnd):
        #vcf_out.write(rec)
        if ct>84:break

        if len(rec.alts) > 1:
            print("Warning found multiallelic variant")
            contine

        print(str(rec.contig) + '\t', str(rec.pos) + '\t', str(rec.ref) + '\t', str(rec.alts[0]) + '\t', end='')

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
            [vf, abhet, passing, failing, missing, gt_failed, depth_sum, clean_obs] = calcVA(sm_list['dict'], rec.filter)

            if subset == 'ADSPccWGS':
                print("VFLAGS_{}={};".format(subset, vf), end='')
                print("ABHet_{}={};".format(subset, abhet), end='')
            write_subset_stats(subset, rec, vf, abhet, passing, failing, missing, gt_failed, depth_sum, clean_obs )

        print()
        ct += 1

    end = time.time()
    print("{0:.2f}".format(end - start))


if __name__ == "__main__":
    main()
    #cProfile.run('main()', None, 'cumtime')
