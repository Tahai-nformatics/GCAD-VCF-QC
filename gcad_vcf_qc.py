#!/usr/bin/env python3
#-*- coding: utf-8 -*-

from argparse import ArgumentParser
from pysam import VariantFile

import csv
import os
import copy
from subprocess import check_output

from collections import namedtuple, OrderedDict, Counter
import time
import cProfile

from utils.VariantAnnotation.vflags import calcVA, calcVA_multiallelic, calcVA_chrx, read_target_files, read_exon_file
from utils.stats.count_gt import is_good_gt, count_gt_multiallelic, count_gt_chrx
from utils.stats.statistical import calc_ExcessHet, calc_ExcessHet_multiallelic, calc_pHWE

import utils.SampleAnnotation.sample_annotation as mi
import config as cfg

import warnings
warnings.simplefilter('always')


def extract_subsets(fam):
    """
    extract_subsets - Creates an OrderedDict() where Subset groupings are the key,
                      values are a set of SampID within it
    @return Samples Dict per Subset
    """

    delimiter = '\t'
    # Check number of columns
    with open(fam, 'r') as fam_file:
      first_line = fam_file.readline()

    ncol = first_line.count(delimiter) + 1
    if ncol == 1:
      delimiter = ','
      ncol = first_line.count(delimiter) + 1

    #
    if ncol == 15:
       SampleFamDetail = namedtuple('SampleFamDetail',['SampID','FID','SubjID','FA','MO',
                                                    'SEX','AFF','AD','AGE','ADSPWGS',
                                                    'Subset','Subgroup','Race_Ethnicity','SeqCtr','ExcludeFromZHet'])
    elif ncol == 17:
       SampleFamDetail = namedtuple('SampleFamDetail',['SampID','FID','SubjID','FA','MO',
                                                    'SEX','AFF','AD','AGE','ADSPWGS',
                                                    'Subset','Subgroup','Race_Ethnicity','SeqCtr','ExcludeFromZHet',
                                                    'TargetFile', 'TargetFilePath'])
    else:
       raise TypeError("Wrong number of columns")


    samples = OrderedDict()
    ct = 0
    with open(fam, 'r') as fam_file:
        for sm in map(SampleFamDetail._make, csv.reader(fam_file, delimiter=delimiter)):
          if ncol==15:
            if sm.Subset in samples:
                samples[sm.Subset].add(sm.SampID)
            else:
                samples[sm.Subset] = set()
                samples[sm.Subset].add(sm.SampID)
          elif ncol==17:
            combined = sm.Subset + '-' + sm.TargetFile

            if combined in samples:
                samples[combined].add(sm.SampID)
            else:
                samples[combined] = set()
                samples[combined].add(sm.SampID)
          ct += 1


    return samples, ct

def extract_subsets_chrx(fam):
    """
    extract_subsets - Creates an OrderedDict() where Subset groupings are the key,
                      values are a set of SampID within it
    @return Samples Dict per Subset
    """

    delimiter = '\t'
    # Check number of columns
    with open(fam, 'r') as fam_file:
      first_line = fam_file.readline()

    ncol = first_line.count(delimiter) + 1
    if ncol == 1:
      delimiter = ','
      ncol = first_line.count(delimiter) + 1

    #
    if ncol == 15:
       SampleFamDetail = namedtuple('SampleFamDetail',['SampID','FID','SubjID','FA','MO',
                                                    'SEX','AFF','AD','AGE','ADSPWGS',
                                                    'Subset','Subgroup','Race_Ethnicity','SeqCtr','ExcludeFromZHet'])
    elif ncol == 17:
       SampleFamDetail = namedtuple('SampleFamDetail',['SampID','FID','SubjID','FA','MO',
                                                    'SEX','AFF','AD','AGE','ADSPWGS',
                                                    'Subset','Subgroup','Race_Ethnicity','SeqCtr','ExcludeFromZHet',
                                                    'TargetFile', 'TargetFilePath'])
    else:
       raise TypeError("Wrong number of columns")

    samples = OrderedDict()
    male_samples = OrderedDict()
    female_samples = OrderedDict()
    ct = 0
    with open(fam, 'r') as fam_file:
        for sm in map(SampleFamDetail._make, csv.reader(fam_file, delimiter=delimiter)):
          if ncol==15:
              if  sm.Subset in male_samples:
                  if sm.SEX == "0":
                    male_samples[sm.Subset].add(sm.SampID)
                  elif sm.SEX == "1":
                    female_samples[sm.Subset].add(sm.SampID)
              else:
                  male_samples[sm.Subset] = set()
                  female_samples[sm.Subset] = set()
                  if sm.SEX == "0":
                    male_samples[sm.Subset].add(sm.SampID)
                  elif sm.SEX == "1":
                    female_samples[sm.Subset].add(sm.SampID)
    
    return male_samples, female_samples



def write_subset_stats_multiallelic(prefix, rec, subset,maf, vf, passing_d,failing_d, missing,gt_failed,clean_passing_d,sum_clean,depth_sum,ab_het,mend_pairs,mend_errors,scores,VTYPE):

    """
         passing_d = {'obs_homo1':{},'obs_het':{},'obs_homo2':{}}
    """
    outfile = '{}.{}.tsv'.format(prefix, subset)
    newfile = not os.path.exists(outfile)
    rec_details = {'filter': rec.filter, 'ref': rec.ref, 'alt': rec.alts,
                 'chr': rec.contig, 'pos': rec.pos}
    total_genotypes = sum_clean + gt_failed


    #MeanDepth    
    mean_depth = depth_sum / total_genotypes if total_genotypes else 0
    #CallRate
    callrate = 1 - (missing + gt_failed) / (missing + gt_failed + sum_clean)

    with open(outfile, 'a') as csvfile:
        fieldnames = ['CHR','POS',
                      'PASS_Homoz_Ref','PASS_Het', "PASS_Homoz_Alt",
                      'FAIL_Homoz_Ref', 'FAIL_Het', 'FAIL_Homoz_Alt',
                      'MISSING', 'GT_FAILED',
                      'CLEAN_Homoz_Ref', 'CLEAN_Het','CLEAN_Homoz_Alt',
                      'MONO','CALLRATE','CALLBAD','GATKPass',
                      'Mendelian_Inconsistency','Mend_pairs','propMI','AF','MEAN_DEPTH', 'HI_DEPTH', 'ABHET',
                      'VFLAGS', 'rsID', 'RefAllele', 'AltAlleles',
                      'QUAL','FILTER','VTYPE',
                     ]
        fieldnames.extend(scores.keys())
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames , delimiter='\t', lineterminator='\n')

        if newfile:
            writer.writeheader()
        qual = "{0:.2f}".format(rec.qual) if rec.qual is not None else 'NA'
        row = {'CHR': rec.contig,
               'POS': rec.pos,
            'PASS_Homoz_Ref':list(passing_d['obs_homo1'].values())[0],
            "PASS_Het":",".join(str(x) for x in passing_d['obs_het'].values()),
            "PASS_Homoz_Alt":",".join(str(x) for x in passing_d['obs_homo2'].values()),
            'FAIL_Homoz_Ref':",".join(str(x) for x in failing_d['obs_homo1'].values()),
            'FAIL_Het':",".join(str(x) for x in failing_d['obs_het'].values()),
            'FAIL_Homoz_Alt':",".join(str(x) for x in failing_d['obs_homo2'].values()),
            'MISSING': missing,
            'GT_FAILED':gt_failed,
            'CLEAN_Homoz_Ref': list(clean_passing_d['obs_homo1'].values())[0],
            'CLEAN_Het':",".join(str(x) for x in clean_passing_d['obs_het'].values()),
            'CLEAN_Homoz_Alt': ",".join(str(x) for x in clean_passing_d['obs_homo2'].values()),
            'MONO': int(3 in vf),
            'CALLRATE':"{0:.5f}".format(callrate), 'CALLBAD':int(callrate < (1 - cfg.miss_rate)),
            'GATKPass': int(1 not in vf),
            'Mendelian_Inconsistency':mend_errors, 'Mend_pairs':mend_pairs,'propMI': "{0:.6f}".format(mend_errors / mend_pairs if mend_pairs >0 else -1),
            'AF': ",".join(str(x) for x in maf),
            'MEAN_DEPTH':"{0:.5f}".format(mean_depth), 'HI_DEPTH':int(mean_depth > cfg.max_dp),
            'ABHET':",".join(str(x) for x in ab_het),
            'VFLAGS': ",".join(str(x) for x in vf),
            'rsID': rec.id if rec.id else '.',
            'RefAllele': rec.ref,
            'AltAlleles': ",".join(rec.alts),
            'QUAL':qual,
            'FILTER':",".join(rec.filter.keys()),
            'VTYPE': VTYPE,
            }
        row.update(scores)
        writer.writerow(row)
    return

def write_subset_stats(prefix, subset, rec, vf, abhet, passing, failing, missing, gt_failed, depth_sum, clean_obs, mend_pairs, mend_errors, scores, vtype, isWES, have_target):
    """
    """
    outfile = '{}.{}.tsv'.format(prefix, subset)
    newfile = not os.path.exists(outfile)

    with open(outfile, 'a') as csvfile:
        fieldnames = ['CHR','POS',
                      'Pass00','Pass01','Pass11',
                      'Fail00','Fail01','Fail11',
                      'Missing','GT_Failed',
                      'Clean00','Clean01','Clean11','Mono','CallRate','CallBad','GATKPass','MAF','AltAF',
                      'MeanDepth','HiDepth','ABHet','Mend_Incon','Mend_pairs','propMI','MultiAllele','FilteredOut',
                      'VFLAGS','rsID','RefAllele','AltAllele','QUAL','FILTER','VTYPE',
                      ]

        # Additional column for WES
        if isWES:
            fieldnames.extend(['InTargetRegion'])

        # Add column names for subgroup scores
        fieldnames.extend(scores.keys())

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames , delimiter='\t', lineterminator='\n')

        if newfile:
            writer.writeheader()

        # field calculations
        sum_clean = sum(clean_obs)

        # callrate
        callrate = 1 - (missing + gt_failed) / (missing + gt_failed + sum_clean)

        # MAF / AltAF; altAF = maf
        maf = 0
        alt_maf = 0
        temp = 2 * sum_clean

        if temp > 0:
            maf = (clean_obs[1] + 2 * clean_obs[2]) / temp
            alt_maf = maf
            if maf > 0.5:
                maf = 1 - maf

        maf = "{0:.6f}".format(maf)

        # MeanDepth
        total_genotypes = sum_clean + gt_failed
        mean_depth = depth_sum / total_genotypes if total_genotypes else 0
        qual = "{0:.2f}".format(rec.qual) if rec.qual is not None else 'NA'


        vf.sort()

        row = {'CHR': rec.contig, 'POS': rec.pos,
            'Pass00':passing[0],'Pass01':passing[1],'Pass11':passing[2],
            'Fail00':failing[0],'Fail01':failing[1],'Fail11':failing[2],
            'Missing':missing, 'GT_Failed':gt_failed,
            'Clean00':clean_obs[0],'Clean01':clean_obs[1],'Clean11':clean_obs[2],
            'Mono':int(3 in vf),
            'CallRate':"{0:.6f}".format(callrate),
            'CallBad': int(callrate < (1 - cfg.miss_rate)),
            'GATKPass':int(1 not in vf),
            'MAF':maf,'AltAF':"{0:.6f}".format(alt_maf),
            'MeanDepth':"{0:.6f}".format(mean_depth),'HiDepth':int(mean_depth > cfg.max_dp),
            'ABHet':abhet,
            'Mend_Incon':mend_errors, 'Mend_pairs':mend_pairs, 'propMI': "{0:.6f}".format(mend_errors / mend_pairs if mend_pairs >0 else 0),
            'MultiAllele':0,'FilteredOut':int(0 not in vf),
            'VFLAGS':",".join(map(str,vf)),
            'rsID':rec.id if rec.id else '.',
            'RefAllele':rec.ref,'AltAllele':",".join(map(str,rec.alts)),
            'QUAL':qual,
            'FILTER':",".join(rec.filter.keys()),
            'VTYPE': vtype
            }

        if isWES:
           if have_target:
              row['InTargetRegion'] = int(11 not in vf)
           else:
              row['InTargetRegion'] = '.'

        row.update(scores)
        writer.writerow(row)

    return

def write_subset_stats_chrx(prefix, rec, subset,vf,passing_d_male,passing_d_female,failing_d_male,failing_d_female,missing,gt_failed,clean_d,sum_clean,maf,depth_sum,ab_het,mend_pairs,mend_errors,scores, vtype):

    """
         passing_d = {GT_type: {GT:count}, GT_type: {GT:count}, GT_type: {GT:count}}

         i.e : passing_d = {'obs_homo1': {((0, 0), (0, 0)): 6309}, 'obs_het': {((0, 1), (1, 0)): 0}, 'obs_homo2': {((1, 1), (1, 1)): 0}}

    """
    rec_details = {'filter': rec.filter, 'ref': rec.ref, 'alt': rec.alts,
                 'chr': rec.contig, 'pos': rec.pos}
    outfile = '{}.{}.tsv'.format(prefix, subset)
    newfile = not os.path.exists(outfile)
    callrate = 1 - (missing + gt_failed) / (missing + gt_failed + sum_clean)

    #MeanDepth
    total_genotypes = sum_clean + gt_failed
    mean_depth = depth_sum / total_genotypes if total_genotypes else 0
    with open(outfile, 'a') as csvfile:
        fieldnames = ['CHR','POS',
                      'PASS_Homoz_Ref','PASS_Het', "PASS_Homoz_Alt",
                      'FAIL_Homoz_Ref', 'FAIL_Het', 'FAIL_Homoz_Alt',
                      'MISSING', 'GT_FAILED',
                      'CLEAN_Homoz_Ref', 'CLEAN_Het','CLEAN_Homoz_Alt',
                      'MONO','CALLRATE','CALLBAD','GATKPass',
                      'AF','MEAN_DEPTH', 'HI_DEPTH', 'ABHET',
                      'Mendelian_Inconsistency','Mend_pairs','propMI',
                      'VFLAGS', 'rsID', 'RefAllele', 'AltAlleles',
                      'QUAL','FILTER','VTYPE','MaleHet',
                     ]
        fieldnames.extend(scores.keys())
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames , delimiter='\t', lineterminator='\n')

        Pass_Het = ",".join(str(x) for x in passing_d_male['obs_het'].values())


        if newfile:
            writer.writeheader()
        qual = "{0:.2f}".format(rec.qual) if rec.qual is not None else 'NA'
        row = {'CHR': rec.contig,
               'POS': rec.pos,
            'PASS_Homoz_Ref':str(list(passing_d_male['obs_homo1'].values())[0])+";"+str(list(passing_d_female['obs_homo1'].values())[0]),
            "PASS_Het":"0"+";"+",".join(str(x) for x in passing_d_female['obs_het'].values()),
            "PASS_Homoz_Alt":",".join(str(x) for x in passing_d_male['obs_homo2'].values())+";"+",".join(str(x) for x in passing_d_female['obs_homo2'].values()),
            'FAIL_Homoz_Ref':",".join(str(x) for x in failing_d_male['obs_homo1'].values())+";"+",".join(str(x) for x in failing_d_female['obs_homo1'].values()),
            'FAIL_Het':",".join(str(x) for x in failing_d_male['obs_het'].values())+";"+",".join(str(x) for x in failing_d_female['obs_het'].values()),
            'FAIL_Homoz_Alt':",".join(str(x) for x in failing_d_male['obs_homo2'].values())+";"+",".join(str(x) for x in failing_d_female['obs_homo2'].values()),
            'MISSING': missing,
            'GT_FAILED':gt_failed,
            'CLEAN_Homoz_Ref': ",".join(str(x) for x in clean_d['male']['obs_homo1'].values()) + ";" + ",".join(str(x) for x in clean_d['female']['obs_homo1'].values()),
            'CLEAN_Het' : str(0) + ";" + ",".join(str(x) for x in clean_d['female']['obs_het'].values()),
            'CLEAN_Homoz_Alt': ",".join(str(x) for x in clean_d['male']['obs_homo2'].values()) + ";" + ",".join(str(x) for x in clean_d['female']['obs_homo2'].values()),
            'MONO': int(3 in vf),
            'CALLRATE':"{0:.5f}".format(callrate), 'CALLBAD':int(callrate < (1 - cfg.miss_rate)),
            'GATKPass': int(1 not in vf),
            'AF': ",".join(str(x) for x in maf),
            'MEAN_DEPTH':"{0:.5f}".format(mean_depth), 'HI_DEPTH':int(mean_depth > cfg.max_dp),
            'ABHET':",".join(str(x) for x in ab_het),
            'Mendelian_Inconsistency':mend_errors, 'Mend_pairs':mend_pairs,'propMI': "{0:.6f}".format(mend_errors / mend_pairs if mend_pairs >0 else -1), 
            'VFLAGS': ",".join(str(x) for x in vf),
            'rsID': rec.id if rec.id else '.',
            'RefAllele': rec.ref,
            'AltAlleles': ",".join(rec.alts),
            'QUAL':qual,
            'FILTER':",".join(rec.filter.keys()),
            'VTYPE': vtype,
            'MaleHet':sum(clean_d['male']['obs_het'].values())
            #'MaleHet': ",".join(str(x) for x in clean_d['male']['obs_het'].values())
            }
        row.update(scores)
        writer.writerow(row)


def write_mendelian_errors(prefix, rec, fam_info, genos ): # mmmm, genos
    """
    """
    outfile = "{}{}".format(prefix, '.tsv')
    newfile = not os.path.exists(outfile)

    with open(outfile, 'a') as csvfile:
        fieldnames = ['FID','PID','CID','CHR','POS','REF','ALT','P1GT','P2GT','CGT',]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames , delimiter='\t', lineterminator='\n')

        if newfile:
            writer.writeheader()

        writer.writerow({'FID': fam_info[0], 'PID':fam_info[1], 'CID': fam_info[2],
                         'CHR': rec.contig, 'POS': rec.pos,
                         'REF': rec.ref, 'ALT':rec.alts[0],
                         'P1GT': "/".join(map(str,genos[0])),
                         'P2GT': "/".join(map(str,genos[1])),
                         'CGT': "/".join(map(str,genos[2]))
                         })


def write_indiv_summary_multiallelic(prefix, isWES):
    """
    """
    outfile = '{}.tsv'.format(prefix)

    with open(outfile, 'w') as csvfile:
        fieldnames = ['SampleID','SEX', 'Pass', 'Fail','Missing', 'Set_Missing','HetHom', 'IndDepthSum', 'IndMeanDepth']

        if isWES:
           fieldnames.extend(['Ti_WES','Tv_WES','TiTvRatio_WES'])

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames , delimiter='\t', lineterminator='\n')

        writer.writeheader()

        abc_order = OrderedDict(sorted(mi.sa.sa_collection.items()))
        for indiv, val in abc_order.items():
            het_hom = val.tallySA['passing_obs_het']/val.tallySA['passing_obs_homo2'] if val.tallySA['passing_obs_homo2'] else 0
            # mean_depth
            good_gt = val.tallySA[(0,0)] + val.tallySA['passing_obs_het'] + val.tallySA['passing_obs_homo2']
            all_gt = good_gt + val.tallySA[-9]
            mean_depth = val.dp_total / all_gt if all_gt else 0
            row = {'SampleID': indiv, 'SEX': val.details_dict.SEX,
                    'Pass': ",".join([str(val.tallySA['passing_obs_homo1']),str(val.tallySA['passing_obs_het']),str(val.tallySA['passing_obs_homo2'])]),
                    'Fail':  ",".join([str(val.tallySA['failing_obs_homo1']),str(val.tallySA['failing_obs_het']),str(val.tallySA['failing_obs_homo2'])]),
                    'Missing': val.tallySA[(None,None)],
                    'Set_Missing': val.tallySA[-9],
                    'HetHom':"{0:.5f}".format(het_hom),
                    'IndDepthSum': val.dp_total,
                    'IndMeanDepth':"{0:.5f}".format(mean_depth),
                    }
            # WES - TiTv
            if isWES:
               ti_tv_wes = val.tallySA['ti_wes'] if val.tallySA['tv_wes'] == 0 else val.tallySA['ti_wes'] / val.tallySA['tv_wes']
               row['Ti_WES'] = val.tallySA['ti_wes']
               row['Tv_WES'] = val.tallySA['tv_wes']
               row['TiTvRatio_WES'] = "{0:.5f}".format(ti_tv_wes)
            writer.writerow(row)
    return


def write_indiv_summary(prefix, isWES):
    """
    """
    outfile = '{}.tsv'.format(prefix)

    with open(outfile, 'w') as csvfile:
        fieldnames = ['SampleID','SEX',
                      'total_nRR','total_nRA','total_nAA','Missing','Set_Missing',
                      'Singleton','Private_Doubleton','Doubleton','HetHom',
                      'Ti','Tv','TiTvRatio','IndDepthSum','IndMeanDepth',
                      '1P_MI','2P_MI','MI_pairs','Non_Missing_Indels',]

        if isWES:
           fieldnames.extend(['Ti_WES','Tv_WES','TiTvRatio_WES'])

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames , delimiter='\t', lineterminator='\n')

        writer.writeheader()

        abc_order = OrderedDict(sorted(mi.sa.sa_collection.items()))
        for indiv, val in abc_order.items():
            # ti/tv
            ti_tv = val.tallySA['ti'] / val.tallySA['tv'] if val.tallySA['tv'] else 0

            # HetHom is the ratio of hets to homozygous alt SNVs
            het_hom = val.tallySA[(0,1)] / val.tallySA[(1,1)] if val.tallySA[(1,1)] else 0

            # mean_depth
            good_het_gt = val.tallySA[(0,1)] + val.tallySA[(1,0)]
            good_gt = val.tallySA[(0,0)] + good_het_gt + val.tallySA[(1,1)]
            all_gt = good_gt + val.tallySA[-9]
            mean_depth = val.dp_total / all_gt if all_gt else 0
            
            row = {'SampleID': indiv, 'SEX': val.details_dict.SEX,
                            'total_nRR': val.tallySA[(0,0)],'total_nRA': good_het_gt,'total_nAA': val.tallySA[(1,1)],
                            'Missing': val.tallySA[(None,None)],'Set_Missing': val.tallySA[-9],
                            'Singleton': val.tallySA['singleton'],
                            'Private_Doubleton': val.tallySA['p_dblton'],
                            'Doubleton': val.tallySA['doubleton'],
                            'HetHom':"{0:.5f}".format(het_hom),
                            'Ti': val.tallySA['ti'], 'Tv': val.tallySA['tv'], 'TiTvRatio':"{0:.5f}".format(ti_tv),
                            'IndDepthSum': val.dp_total,
                            'IndMeanDepth':"{0:.5f}".format(mean_depth),
                            '1P_MI': val.tallySA['vp1'],'2P_MI':val.tallySA['vp2'],'MI_pairs':val.tallySA['mend_pair'],
                            'Non_Missing_Indels': val.tallySA['non_missing_indel']
                            }

            # WES - TiTv
            if isWES:
               ti_tv_wes = val.tallySA['ti_wes'] if val.tallySA['tv_wes'] == 0 else val.tallySA['ti_wes'] / val.tallySA['tv_wes']
               row['Ti_WES'] = val.tallySA['ti_wes']
               row['Tv_WES'] = val.tallySA['tv_wes']
               row['TiTvRatio_WES'] = "{0:.5f}".format(ti_tv_wes)
            writer.writerow(row)
    return


def write_indiv_summary_chrx(prefix, isWES, is_multiallelic):
    """
    """
    outfile = '{}.tsv'.format(prefix)
    with open(outfile, 'w') as csvfile:
        
        if is_multiallelic:
            fieldnames = ['SampleID','SEX', 'Pass', 'Fail','Missing', 'Set_Missing',
                      'HetHom','IndDepthSum','IndMeanDepth']
        else:
            fieldnames = ['SampleID','SEX', 'Pass', 'Fail','Missing', 'Set_Missing',
                      'Singleton','Private_Doubleton','Doubleton','HetHom',
                      'Ti','Tv','TiTvRatio','IndDepthSum', 'IndMeanDepth',]

        if isWES: # Deprecated
           fieldnames.extend(['Ti_WES','Tv_WES','TiTvRatio_WES'])

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames , delimiter='\t', lineterminator='\n')
        writer.writeheader()
        abc_order = OrderedDict(sorted(mi.sa.sa_collection.items()))
        for indiv, val in abc_order.items():
            if not is_multiallelic:
                ti_tv = val.tallySA['ti'] / val.tallySA['tv'] if val.tallySA['tv'] else 0

            het_hom = (val.tallySA['passing_obs_het'] /val.tallySA['passing_obs_homo2'] if val.tallySA['passing_obs_homo2'] else 0) if val.details_dict.SEX == "1" else 0

            # mean_depth
            good_gt = val.tallySA[(0,0)] + val.tallySA['passing_obs_homo2'] + val.tallySA['passing_obs_het']
            all_gt = good_gt + val.tallySA[-9]
            mean_depth = val.dp_total / all_gt if all_gt else 0
            #Multiallelic ChrX Rows
            if is_multiallelic:
                row = {'SampleID': indiv, 'SEX': val.details_dict.SEX,
                    'Pass': ",".join([str(val.tallySA['passing_obs_homo1']),str(val.tallySA['passing_obs_het']) if val.details_dict.SEX == "1" else "0",str(val.tallySA['passing_obs_homo2'])]),
                    'Fail':  ",".join([str(val.tallySA['failing_obs_homo1']),str(val.tallySA['failing_obs_het']),str(val.tallySA['failing_obs_homo2'])]),
                    'Missing': val.tallySA[(None,None)],
                    'Set_Missing': val.tallySA[-9],
                    'HetHom':"{0:.5f}".format(het_hom),
                    'IndDepthSum': val.dp_total,
                    'IndMeanDepth':"{0:.5f}".format(mean_depth),
                    }
            else: #Biallelic Chrx Rows
                row = {'SampleID': indiv, 'SEX': val.details_dict.SEX,
                    'Pass': ",".join([str(val.tallySA['passing_obs_homo1']),str(val.tallySA['passing_obs_het']) if val.details_dict.SEX == "1" else "0",str(val.tallySA['passing_obs_homo2'])]),
                    'Fail':  ",".join([str(val.tallySA['failing_obs_homo1']),str(val.tallySA['failing_obs_het']),str(val.tallySA['failing_obs_homo2'])]),
                    'Missing': val.tallySA[(None,None)],
                    'Set_Missing': val.tallySA[-9],
                    'Singleton': val.tallySA['singleton'],
                    'Private_Doubleton': val.tallySA['p_dblton'],
                    'Doubleton': val.tallySA['doubleton'],
                    'HetHom':"{0:.5f}".format(het_hom),
                    'Ti': val.tallySA['ti'], 
                    'Tv': val.tallySA['tv'], 
                    'TiTvRatio':"{0:.5f}".format(ti_tv),
                    'IndDepthSum': val.dp_total,
                    'IndMeanDepth':"{0:.5f}".format(mean_depth),
                    }

            # WES - TiTv
            if isWES:
               ti_tv_wes = val.tallySA['ti_wes'] if val.tallySA['tv_wes'] == 0 else val.tallySA['ti_wes'] / val.tallySA['tv_wes']
               row['Ti_WES'] = val.tallySA['ti_wes']
               row['Tv_WES'] = val.tallySA['tv_wes']
               row['TiTvRatio_WES'] = "{0:.5f}".format(ti_tv_wes)
            writer.writerow(row)
    return


def delete_previous_outputs(out_dir, prefix, subsets):
    """
    """
    mi_file = out_dir + 'summary.mi.tsv'
    for out_file in [mi_file] + [ '{}{}.{}.tsv'.format(out_dir, prefix, x) for x in subsets]:
        if os.path.exists(out_file):
            os.remove(out_file)
    return



def vcf_output_create_biallelic(rec, subset, clean_obs, vf, abhet, vtype, vcf_out):
    sum_clean = sum(clean_obs)
    maf = 0
    alt_maf = 0
    temp = 2 * sum_clean
    alt_allele_counts = 2*clean_obs[2] + clean_obs[1]

    if temp > 0:
        maf = (clean_obs[1] + 2 * clean_obs[2]) / temp
        alt_maf = maf
        if maf > 0.5:
            maf = 1 - maf

        maf = "{0:.6f}".format(maf)

    #Append Allele Number to INFO field
    rec.info["AN"] = sum_clean
    

    #Append Allele Counts to INFO field
    print(rec.info['AC'])
    rec.info['AC'] = alt_allele_counts
    #print(rec.info['AC'])

    #Append Allele Frequency to INFO field
    #rec.info['AF'] = tuple(float(maf[key]) for key in alt_allele_counts.keys())
    rec.info['AF'] = float(alt_maf)

    #Append VFLAGS to INFO field
    rec.info[ "VFLAGS_" + subset ] = vf
    #Append subset ABHet to INFO field
    rec.info[ "ABHet_" + subset ] = [str(num) for num in abhet]

    #Append VariantType to INFO field 
    rec.info[ "VariantType" ] = vtype

    #Write to VCF File
    vcf_out.write(rec)




def vcf_output_create_multiallelic(rec, subset, clean_d, maf, vf, ab_het, vtype, vcf_out):
    total_sum = sum(sum(Genotype_type.values()) for Genotype_type in clean_d.values()) *2
    alt_allele_counts = {i: 0 for i in range(1,len(rec.alts)+1)}

    #Loop through Passing dictionary and add allele counts to alt_allele_counts (For AC)
    for allele in alt_allele_counts:
        for GT_type in clean_d.keys():
            for key,values in clean_d[GT_type].items():
                if any(allele in subkey for subkey in key):
                    if GT_type == "obs_het":
                        alt_allele_counts[allele] += clean_d[GT_type][key]
                    else:
                        alt_allele_counts[allele] += (clean_d[GT_type][key] *2)

    #Append Allele Number to INFO field
    rec.info["AN"] = total_sum
    
    #Append Allele Counts to INFO field
    rec.info['AC'] = tuple(alt_allele_counts[key] for key in alt_allele_counts.keys())
    
    #Append Allele Frequency to INFO field
    rec.info['AF'] = tuple(float(maf[key]) for key in alt_allele_counts.keys())

    #Append VFLAGS to INFO field
    rec.info[ "VFLAGS_" + subset ] = vf

    #Append subset ABHet to INFO field
    rec.info[ "ABHet_" + subset ] = [str(num) for num in ab_het]

    #Append VariantType to INFO field 
    rec.info[ "VariantType" ] = vtype

    #Write to VCF File
    vcf_out.write(rec)

def vcf_output_create_chrX(rec, subset, passing_d, maf, vf, ab_het, vtype, vcf_out, chrx_is_multiallelic):
    alt_allele_counts = {i: 0 for i in range(1,len(rec.alts)+1)}
    total_sum_male = 0
    total_sum_female = 0
    for sex in passing_d.keys():
        if sex == 'male':
            for key, values in passing_d['male'].items():
                if key == 'obs_homo1' or key == 'obs_homo2':
                    total_sum_male += sum(values.values())
        else: #Female
            for key,values in passing_d['female'].items():
                total_sum_female += sum(values.values())
    total_sum = (total_sum_female *2) + ( total_sum_male )

    #Loop through Passing dictionary and add allele counts to alt_allele_counts (For AC)
    for allele in alt_allele_counts:
        for sex in passing_d.keys():
            for GT_type in passing_d[sex].keys():
                if sex == 'male' and GT_type != 'obs_homo2': 
                    continue
                for key,values in passing_d[sex][GT_type].items():
                    if any(allele in subkey for subkey in key):
                        if GT_type == "obs_het":
                            alt_allele_counts[allele] += (passing_d['female'][GT_type][key])
                        elif sex == "male": #Can only be obs_homo2
                            alt_allele_counts[allele] += (passing_d['male'][GT_type][key])
                        else: #Female obs_homo2
                            alt_allele_counts[allele] += (passing_d['female'][GT_type][key] *2)
    
    #Append Allele Number to INFO field
    rec.info["AN"] = total_sum
    
    #Append Allele Counts to INFO field
    rec.info['AC'] = tuple(alt_allele_counts[key] for key in alt_allele_counts.keys())
    
    #Append Allele Frequency to INFO field
    rec.info['AF'] = tuple(float(maf[key]) for key in alt_allele_counts.keys())
    
    #Append VFLAGS to INFO field
    rec.info[ "VFLAGS_" + subset ] = vf

    #Append subset ABHet to INFO field
    #abhet_vcf = [str(item) for item in ab_het]
    rec.info[ "ABHet_" + subset ] = [str(num) for num in ab_het]

    #Append VariantType to INFO field
    rec.info[ "VariantType" ] = vtype

    #Write to VCF File
    vcf_out.write(rec)


def main():
    argparser = ArgumentParser()
    grp_file_paths = argparser.add_argument_group(title='File paths')
    grp_file_paths.add_argument('--vcf', type=str, help='input VCF file', required=True)
    grp_file_paths.add_argument('--fam', type=str, help='fam file (with headers)', required=True)
    grp_file_paths.add_argument('--out_dir', type=str, help='filepath for output file', default='./', required=False)

    grp_overrides = argparser.add_argument_group(title='BCF Overrides')
    grp_overrides.add_argument('--region', type=str, help='restrict to VCF region e.g. chr3:100-200', default=None, required=False)

    fun_overrides = argparser.add_argument_group(title='Functional Overrides')
    fun_overrides.add_argument('--no_output_vcf', default=False, action="store_true", help='use to skip output VCF', required=False)

    grp_settings = argparser.add_argument_group(title='Threshold Settings')
    grp_settings.add_argument('--min_dp', type=int, help='minimum depth (DP)', default=10, required=False)
    grp_settings.add_argument('--min_gq', type=int, help='minimum genotype quality (GQ)', default=20, required=False)
    grp_settings.add_argument('--is_fam', type=int, help='is this family data, 1=yes, 0=no', default=1, required=False)
    grp_settings.add_argument('--min_tranche', type=int, help='Minimum Tranche Score (Filter from VQSR)', default=99.7, required=False)
    grp_settings.add_argument('--miss_rate', type=int, help='max missingness threshold', default=0.2, required=False)
    grp_settings.add_argument('--max_dp', type=int, help='max DP threshold', default=500, required=False)
    grp_settings.add_argument('--hetz_lim1', type=int, help='', default=99999, required=False)
    grp_settings.add_argument('--hetz_lim2', type=int, help='', default=99999, required=False)
    grp_settings.add_argument('--hwe_pval', type=int, help='', default=5e-06, required=False)
    grp_settings.add_argument('--hwe_maf', type=int, help='MAF threshold', default=0.01, required=False)
    grp_settings.add_argument('--is_multiallelic', help='is this vcf multiallelic', action='store_true', required=False)
    grp_settings.add_argument('--is_chrx', help='is this vcf chrx', action='store_true', required=False)

    wes_settings = argparser.add_argument_group(title='WES Settings')
    wes_settings.add_argument('--wes_target', type=str, help='WES Target BED file, e.g. filename:subset ', nargs='*', required=False)
    wes_settings.add_argument('--flank_size', type=int, help='size of target interval expansion in bp', default=7, required=False)
    wes_settings.add_argument('--exon_file',  type=str, help='exon file for TiTv_WES', default=None, required=False)
    wes_settings.add_argument('--chr',        type=str, help='chromosome loading hint', default=None, required=False)

    args, xtra = argparser.parse_known_args()
    print(args)
    print(xtra)

    start = time.time()
    rChr = None
    rStart = None
    rEnd = None
    regionStr = ''

    if args.region:
        rChr = args.region.split(':')[0]
        rStart = int(args.region.split(':')[1].split('-')[0])
        rEnd = int(args.region.split(':')[1].split('-')[1])
        regionStr = ".{}:{}-{}".format(rChr, rStart, rEnd)
        rStart -= 1
        if rStart < 0: rStart = 0
    if rChr is None and args.chr is not None:
       rChr = args.chr
    # Setup globals
    cfg.MINDP = args.min_dp
    cfg.MINGQ = args.min_gq
    cfg.minTranche = args.min_tranche
    cfg.miss_rate = args.miss_rate
    cfg.max_dp = args.max_dp
    cfg.isFam = args.is_fam
    cfg.hetz_lim1 = args.hetz_lim1
    cfg.hetz_lim2 = args.hetz_lim2
    cfg.hwe_pval = args.hwe_pval
    cfg.hwe_maf = args.hwe_maf
    cfg.flank_size = args.flank_size

    # setup output_dir
    if not args.out_dir.endswith("/"): args.out_dir += "/"

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir, exist_ok = True)

    # Create global sample data in mi
    mi.createSampleAnnotation(args.fam)
    samplesDict, famCt = extract_subsets(args.fam)  # returns dict

    print("[FAM] Found {} subsets: {}; totaling {} sampIDs".format(len(samplesDict.keys()), list(samplesDict.keys()), famCt))
    print("[FAM] {}".format([  "{}:{}".format(k, len(samplesDict[k]))  for k in samplesDict.keys()]))
    print("[FAM] {} unique subgroups:{}".format(len(mi.sa.subgroups),sorted( mi.sa.subgroups)))
    print("[FAM] {} ".format([  "{}:{}".format(k, v)  for k,v in mi.sa.subsets.items() ]))

    # INPUT VCF
    vcf_in = VariantFile(args.vcf)
    set_size_in = len(vcf_in.header.samples)

    # Save groups of samples as an intersecting set
    set_size = 0
    for k,v in samplesDict.items():
        samplesDict[k] = {'set': set(vcf_in.header.samples) & v,'dict':dict() }
        set_size += len(samplesDict[k]['set'])

    # Read only a subset of samples
    if set_size != set_size_in:
        vcf_in.subset_samples(mi.sa.id_list)

    # Process WES TARGET BED(s)
    isWES = 0
    total_targets = 0
    if args.wes_target and len(args.wes_target) > 0:
        isWES = 1
        print("[WES] Reading-in target interval files for subsets: {}".format(args.wes_target))
        read_target_files(args.wes_target, rChr)

        # Determine WES target file presence
        for subset in samplesDict.keys():
            have_target = 0
            for tgt in args.wes_target:
               if subset in tgt:
                   have_target = 1
                   total_targets += 1

            samplesDict[subset]['have_target'] = have_target

    if mi.sa.get_targets():
        isWES = 1
#       print("[WES] Reading-in target interval files for samples")
        read_target_files(mi.sa.get_targets(), rChr)

        # Determine WES target file presence
        for subset in samplesDict.keys():
            have_target = 0
            for tgt in mi.sa.get_targets():
               if subset in tgt:
                   have_target = 1
                   total_targets += 1

            samplesDict[subset]['have_target'] = have_target

    # Process exon file
    if args.exon_file:
       print("[WES] Reading Exons BED file")
       read_exon_file(args.exon_file, rChr)

    # organize new vcf_out header
    vcf_out_hdr = vcf_in.header


    for k in samplesDict.keys():
        vcf_out_hdr.add_meta('INFO', items=[('ID', 'VFLAGS_' + k), ('Number','.'), ('Type', 'Integer'), ('Description','Pipeline-specific QC variant flags')])
        if args.is_multiallelic or args.is_chrx:
            vcf_out_hdr.add_meta('INFO', items=[('ID', 'ABHet_' + k), ('Number','.'), ('Type', 'String'), ('Description','Allelic Read Ratio')])
        else: #Biallelic
            vcf_out_hdr.add_meta('INFO', items=[('ID', 'ABHet_' + k), ('Number','1'), ('Type', 'Float'), ('Description','Allelic Read Ratio')])


            
    if isWES:
        vcf_out_hdr.add_meta('INFO', items=[('ID', 'VariantInTargetFraction'), ('Number','.'), ('Type', 'String'), ('Description','Fraction of the variant\'s presence in given target regions')])
        vcf_out_hdr.add_meta('INFO', items=[('ID', 'VariantInTargetRatio'), ('Number',1), ('Type', 'Float'), ('Description','Ratio of the variant\'s presence in given target regions')])

    vcf_out_hdr.add_meta('INFO', items=[('ID', 'VariantType'), ('Number',1), ('Type', 'String'), ('Description','Variant type description')])

    vcf_out_hdr.add_meta('qc_tool', value=os.path.basename(__file__))
    vcf_out_hdr.add_meta('qc_tool-version', value=check_output(["git", "--git-dir", os.path.dirname(__file__) + "/.git", "rev-parse", "--short", "HEAD"]).strip())
    vcf_out_hdr.add_meta('qc_tool-arguments', value="{}".format(args))

    # Output filename for companions
    prefix_companions = args.out_dir + 'summary.snv' + regionStr
    # Output filename for indiv summary
    prefix_indiv = args.out_dir + 'summary.indiv' + regionStr
    # Output filename for MI
    prefix_mi = args.out_dir + 'summary.mi' + regionStr

    print("[IN VCF] contains {} samples".format(set_size_in))

    # Output VCF
    if args.no_output_vcf:
        print("VCF output disabled.")
        print("{}".format(args))
    else:
        baseStr = os.path.basename(args.vcf.replace('.g.vcf','').replace('.vcf','').rpartition('.')[0])
        vcf_out_filename = "{}{}".format(args.out_dir, 'flagged.' + baseStr + regionStr + '.g.vcf.gz')
        vcf_out = VariantFile(vcf_out_filename, 'w', header=vcf_out_hdr, threads=2)

        print("[OUT VCF] will have {} samples from intersecting set".format(set_size))

    delete_previous_outputs(args.out_dir, 'summary.snv' + regionStr, list(samplesDict.keys()))

    ct = 0
    start_p = time.time()
    
    ## Run analysis on Multiallelic chromosome ##
    if args.is_multiallelic:
        if not args.is_chrx:
            for rec in vcf_in.fetch(rChr, rStart, rEnd):
                #Calculate Vtype
                alts_vtype = []
                for alts in rec.alts:
                    ct = 0
                    if alts != "*":
                        if len(rec.ref) != len(alts):
                            vtype='INDEL'
                            alts_vtype.append(vtype)
                            continue
                        for i in range(len(rec.ref)):   # IF length of ref and alt is same...
                            if rec.ref[i] != alts[i]:   #If there is mismatch in BP
                                vtype='SNP'
                                ct +=1           #
                        if ct >1:
                            vtype='INDEL'
                            alts_vtype.append(vtype)
                        elif ct == 1:
                            VTYPE='SNP'
                            alts_vtype.append(vtype)

                if ( 'SNP' in alts_vtype ) and ( 'INDEL' in alts_vtype ):
                    vtype='MULTI_MIX'
                elif 'SNP' in alts_vtype:
                    vtype='MULTI_SNP'
                elif 'INDEL' in alts_vtype:
                    vtype='MULTI_INDEL'
                else:
                    raise "VTYPE error" 

                samplesDict = gather_intersect_fam_vcf_samples(rec.samples, samplesDict)
                for subset, sm_list in samplesDict.items():
                    
                    rec_details = {'filter': rec.filter, 'ref': rec.ref, 'alt': rec.alts,
                                             'chr': rec.contig, 'pos': rec.pos}
                    #Calculate stats
                    [vf,maf,passing_d,failing_d,missing,gt_failed,clean_passing_d,sum_clean,depth_sum,ab_het,subg,subg_c,allele_count_dict,zhet_sample_counts] = calcVA_multiallelic(sm_list['dict'],rec_details,subset)
                    
                    #MI
                    mend_pairs, mend_errors = check_mendelian_errors_multiallelic(prefix_mi, rec)
                    
                    #pHWE per subgroup
                    scores = calculate_subgroup_scores_multiallelic(subset, subg, subg_c,allele_count_dict,zhet_sample_counts)
                    
                    #Companion file
                    write_subset_stats_multiallelic(prefix_companions, rec, subset, maf ,vf,passing_d,failing_d,missing,gt_failed,clean_passing_d,sum_clean,depth_sum,ab_het,mend_pairs,mend_errors,scores,vtype)
                    

                find_s_d_multiallelic(clean_passing_d,rec.samples)
           
                if args.no_output_vcf == False:
                    #Append to INFO field headers and write to VCF file
                    vcf_output_create_multiallelic(rec, subset, clean_passing_d, maf, vf, ab_het, vtype, vcf_out)

                ct += 1
            
            if args.no_output_vcf == False:
                vcf_out.close()

            end = time.time()
            print("total_time:{0:.2f}".format(end - start))
            print("process_time:{0:.2f}".format(end - start_p))
            print("total_processed:{}".format(ct))
            if ct > 0:
                print("rate:{0:.1f}".format(ct/(end - start_p)))
            write_indiv_summary_multiallelic(prefix_indiv, isWES)

        if args.no_output_vcf == False:
            # create index
            time.sleep(1)
            check_output(["tabix", "-f", vcf_out_filename])

    ## Run analysis on ChrX (Biallelic/Multiallelic) chromosome ##
    if args.is_chrx:
        samplesDict_male,samplesDict_female = extract_subsets_chrx(args.fam)

        for k,v in samplesDict_male.items():
            samplesDict_male[k] = {'set': set(vcf_in.header.samples) and v,'dict':dict() }
            set_size += len(samplesDict_male[k]['set'])
            for k,v in samplesDict_female.items():
                samplesDict_female[k] = {'set': set(vcf_in.header.samples) and v,'dict':dict() }
                set_size += len(samplesDict_female[k]['set'])
        for rec in vcf_in.fetch(rChr, rStart, rEnd):
            if len(rec.alts) > 1:
            #Initiate multiallelic variable for write_indiv_summary_chrx and find vtype
                chrx_is_multiallelic = True
                alts_vtype = []
                for alts in rec.alts:
                    ct = 0
                    if alts != "*":
                        if len(rec.ref) != len(alts):
                            vtype='INDEL'
                            alts_vtype.append(vtype)
                            continue
                        for i in range(len(rec.ref)):   # IF length of ref and alt is same...
                            if rec.ref[i] != alts[i]:   #If there is mismatch in BP
                                vtype='SNP'
                                ct +=1           #
                        if ct >1:
                            vtype='INDEL'
                            alts_vtype.append(vtype)
                        elif ct == 1:
                            vtype='SNP'
                            alts_vtype.append(vtype)

                if ( 'SNP' in alts_vtype ) and ( 'INDEL' in alts_vtype ):
                    vtype='MULTI_MIX'
                elif 'SNP' in alts_vtype:
                    vtype='MULTI_SNP'
                elif 'INDEL' in alts_vtype:
                    vtype='MULTI_INDEL'
                else:
                    raise "VTYPE error" 
             
            else: #Initiate biallelic variable for write_indiv_summary_chrx and find vtype
                chrx_is_multiallelic = False
                vtype = "SNV"
                if len(rec.ref) > 1:
                    vtype = "Deletion"
                elif len(rec.alts[0]) > 1:
                    vtype = "Insertion" 

            samplesDict_male = gather_intersect_fam_vcf_samples(rec.samples, samplesDict_male)
            samplesDict_female = gather_intersect_fam_vcf_samples(rec.samples, samplesDict_female)
            for (subset_male, sm_list_male), (subset_female, sm_list_female) in zip(samplesDict_male.items(), samplesDict_female.items()):
                
                rec_details= {'filter': rec.filter, 'ref': rec.ref, 'alt': rec.alts,
                     'chr': rec.contig, 'pos': rec.pos}
                
                #Calculate stats
                [vf,passing_d_male,passing_d_female,failing_d_male,failing_d_female,missing,gt_failed,clean_d,sum_clean,maf,depth_sum, ab_het,subg_male,subg_female,subg_c_male,subg_c_female,zhet_dict,zhet_sample_counts] = calcVA_chrx(sm_list_male['dict'],sm_list_female['dict'],rec_details,subset_male,subset_female)
                #MI 
                mend_pairs, mend_errors = check_mendelian_errors_chrx(prefix_mi, rec)
                
                #pHWE per subgroup
                scores = calculate_subgroup_scores_chrx(rec.alts,subset_male, subg_male,subg_female, subg_c_male,subg_c_female,zhet_dict,zhet_sample_counts)
                
                #Companion file
                write_subset_stats_chrx(prefix_companions, rec, subset_male,vf,passing_d_male,passing_d_female,failing_d_male,failing_d_female,missing,gt_failed,clean_d,sum_clean,maf,depth_sum,ab_het,mend_pairs, mend_errors, scores, vtype) 

            find_s_d_chrx(clean_d, rec.samples)

            if args.no_output_vcf == False:
                #Append to INFO field headers and write to VCF file
                vcf_output_create_chrX(rec, subset_female, clean_d, maf, vf, ab_het, vtype, vcf_out, chrx_is_multiallelic)
                #vcf_out.write(rec)

            ct += 1

        if args.no_output_vcf == False:
            vcf_out.close()

        end = time.time()
        print("total_time:{0:.2f}".format(end - start))
        print("process_time:{0:.2f}".format(end - start_p))
        print("total_processed:{}".format(ct))
        if ct > 0:
            print("rate:{0:.1f}".format(ct/(end - start_p)))
        
        write_indiv_summary_chrx(prefix_indiv, isWES, chrx_is_multiallelic)
    
        if args.no_output_vcf == False:
        # create index
            time.sleep(1)
            print('here')
            check_output(["tabix", "-f", vcf_out_filename])
    

    ## Run analysis on Biallelic Chromosome ##
    elif not args.is_multiallelic and not args.is_chrx: 
        # loop over each variant in VCF
        for rec in vcf_in.fetch(rChr, rStart, rEnd):
            if len(rec.alts) > 1:
            #print("Warning found multiallelic variant")
                continue
        # Variant Type: SNV, MNV, insertion, deletion
            vtype = "SNV"
            if len(rec.ref) > 1:
                vtype = "Deletion"
            elif len(rec.alts[0]) > 1:
                vtype = "Insertion"

            #
            samplesDict = gather_intersect_fam_vcf_samples(rec.samples, samplesDict)
            grp_obs = list()
            vflag_11_ct = 0

            for subset, sm_list in samplesDict.items():
                # calc stats
                [vf, abhet, passing, failing, missing, gt_failed, depth_sum, clean_obs, subg, subg_cntl] = calcVA(
                    sm_list['dict'],
                    {'filter': rec.filter, 'ref': rec.ref, 'alt': rec.alts,
                     'chr': rec.contig, 'pos': rec.pos
                     },
                    subset,
                )

                grp_obs.append(clean_obs)

                # MI
                mend_pairs, mend_errors = check_mendelian_errors(prefix_mi, rec)

                # pHWE per subgroup
                scores = calculate_subgroup_scores(subset, subg, subg_cntl)

                # Companion file
                have_target = 0
                if isWES:
                   have_target = samplesDict[subset]['have_target']
                write_subset_stats(prefix_companions, subset, rec, vf, abhet,
                                   passing, failing, missing, gt_failed, depth_sum, clean_obs,
                                   mend_pairs, mend_errors, scores, vtype,
                                   isWES, have_target
                                  )
      

                #vcf_output_create_biallelic(rec, subset, clean_obs, vf, abhet, vtype, vcf_out)
                
                #Append Allele Number to INFO field
                rec.info["AN"] =  2 * sum(clean_obs)
                
                #Append Allele Counts to INFO field
                rec.info['AC'] = 2*clean_obs[2] + clean_obs[1]
                
                #Append Allele Frequency to INFO field
                rec.info['AF'] = float((clean_obs[1] + (2 * clean_obs[2]))/ (2 * sum(clean_obs)))

                # Append subset VFLAGS to INFO field
                rec.info[ "VFLAGS_" + subset ] = vf

                # Append subset ABHet to INFO field
                rec.info[ "ABHet_" + subset ] = float(abhet) if abhet != 'NA' else None

                # Append VariantType
                rec.info[ "VariantType" ] = vtype

                # count number of vflag(11) for VariantInTargetRatio
                if have_target:
                   #vflag_11_ct += (11 not in vf)
                   if (11 not in vf):
                      #vflag_11_ct += sum( mi.sa.subsets[subset].values())
                      vflag_11_ct += missing + gt_failed + sum(clean_obs)

            # Append VariantInTargetRatio
            if isWES:
               rec.info[ "VariantInTargetFraction" ] = str(vflag_11_ct) + '/' + str(set_size)
               rec.info[ "VariantInTargetRatio" ] = vflag_11_ct / set_size

            # sum obs by column
            total_obs = list(map(sum, zip(*grp_obs)))

            find_s_d(total_obs, rec.samples)

            if args.no_output_vcf == False:
                vcf_out.write(rec)

            #print()
            ct += 1

        if args.no_output_vcf == False: vcf_out.close()

        end = time.time()
        print("total_time:{0:.2f}".format(end - start))
        print("process_time:{0:.2f}".format(end - start_p))
        print("total_processed:{}".format(ct))
        if ct > 0:
            print("rate:{0:.1f}".format(ct/(end - start_p)))
        write_indiv_summary(prefix_indiv, isWES)

        if args.no_output_vcf == False:
            # create index
            time.sleep(1)
            check_output(["tabix", "-f", vcf_out_filename])

def calculate_subgroup_scores_multiallelic(subset, subg, subg_cntl,allele_count_dict,zhet_sample_counts):
    """ calculate_subgroup_scores - generates nClean, Zhet, and pHWE for subgroups
                                    added to TAGs within the INFO field. pHWE-subgroup has
                                    the following criteria, (1) must have N >= 5,
                                    (2) must only use data from controls within the subgroup

        @return scores - dict() of the added calculations
    """
    from itertools import combinations_with_replacement
    scores = OrderedDict()
    
    if mi.sa.get_divide():
      subset = subset.split('-')[0]
    for k in sorted(mi.sa.subsets[subset]):
        zhet_val = allele_count_dict[k]
        zhet_count = zhet_sample_counts[k][0]
        phwe_vals = {k:[0,0,0]}
        nclean_subg = {k:[]}
        nclean_cntl = {k:[]}
        phwe_vals = {k:[0,0,0]}
        total_obs = 0
        alts = len(allele_count_dict[k]) -1

        for gts in subg_cntl[k].keys():
            total_obs += sum([x + y for x,y in zip(subg[k][gts].values(),subg_cntl[k][gts].values())])
            
        for classification in subg_cntl[k]:
            if classification == 'obs_homo1':
                phwe_vals[k][0] = sum(subg_cntl[k][classification].values())
            elif classification == 'obs_het':
                for key,val in subg_cntl[k][classification].items():
                    if not 0 in key[0]:
                        phwe_vals[k][2] += subg_cntl[k][classification][key]
                    else:
                        phwe_vals[k][1] += subg_cntl[k][classification][key]
            else:
                phwe_vals[k][2] += sum(subg_cntl[k][classification].values())
   

        genotypes = combinations_with_replacement(list(n for n in range(alts+1)), 2)
        for allele in genotypes:
            allele_match = (allele[0], allele[1])
            for classification in subg[k]:
                for (key1, value1), (key2,value2) in zip(subg[k][classification].items(), subg_cntl[k][classification].items()):
                    if allele_match in key1:
                        nclean_subg[k].append(subg[k][classification][key1])
                        nclean_cntl[k].append(subg_cntl[k][classification][key2])

        scores['nClean_' + k] = ",".join((str(x) for x in nclean_subg[k])) + ';' +  ",".join((str(x) for x in nclean_cntl[k]))
        scores['Zhet_' + k] = calc_ExcessHet_multiallelic(zhet_val,zhet_count)
        scores['pHWE_' + k] = calc_pHWE(*phwe_vals[k]) if sum(phwe_vals[k]) >= 5 else '.'
        if type(scores['Zhet_' + k]) == float:
            scores['Zhet_' + k] = "{0:.6f}".format(scores['Zhet_' + k])
        if type(scores['pHWE_' + k]) == float:
            scores['pHWE_' + k] = "{0:.6f}".format(scores['pHWE_' + k])
    

    return scores


def calculate_subgroup_scores(subset, subg, subg_cntl):
    """ calculate_subgroup_scores - generates nClean, Zhet, and pHWE for subgroups
                                    added to TAGs within the INFO field. pHWE-subgroup has
                                    the following criteria, (1) must have N >= 5,
                                    (2) must only use data from controls within the subgroup

        @return scores - dict() of the added calculations
    """
    scores = OrderedDict()

    if mi.sa.get_divide():
      subset = subset.split('-')[0]

    for k in sorted(mi.sa.subsets[subset]):
        val = subg[k]
        val_cntl = subg_cntl[k]
        scores['nClean_' + k] =  ",".join(map(str,val)) + ';' + ",".join(map(str,val_cntl))
        scores['Zhet_' + k] = calc_ExcessHet(*val)[0]
        scores['pHWE_' + k] = calc_pHWE(*val_cntl) if sum(val_cntl) >= 5 else '.'
        if type(scores['Zhet_' + k]) == float:
            scores['Zhet_' + k] = "{0:.6f}".format(scores['Zhet_' + k])
        if type(scores['pHWE_' + k]) == float:
            scores['pHWE_' + k] = "{0:.6f}".format(scores['pHWE_' + k])

    return scores

def calculate_subgroup_scores_chrx(alts,subset, subg_male,subg_female, subg_cntl_male,subg_cntl_female,zhet_dict,zhet_sample_counts):
    """ calculate_subgroup_scores - generates nClean, Zhet, and pHWE for subgroups
                                    added to TAGs within the INFO field. pHWE-subgroup has
                                    the following criteria, (1) must have N >= 5,
                                    (2) must only use data from controls within the subgroup

        @return scores - dict() of the added calculations
    """
    from itertools import combinations_with_replacement
    scores = OrderedDict()

    if mi.sa.get_divide():
      subset = subset.split('-')[0]

    for k in sorted(mi.sa.subsets[subset]):
        zhet_val = zhet_dict[k]
        nclean_female_subg = {k:[]}
        nclean_female_subg_cntl = {k:[]}
        phwe_vals = {k:[0,0,0]}
        zhet_count = zhet_sample_counts[k][0]

        genotypes = combinations_with_replacement(list(n for n in range(len(alts)+1)), 2)
        for allele in genotypes:
            allele_match = (allele[0], allele[1])
            for classification in subg_female[k]:
                for (key1, value1), (key2,value2) in zip(subg_female[k][classification].items(), subg_cntl_female[k][classification].items()):
                    if allele_match in key1:
                        nclean_female_subg[k].append(subg_female[k][classification][key1])
                        nclean_female_subg_cntl[k].append(subg_cntl_female[k][classification][key2])


        #Treat multiallelic Heterozygous controls with non 0 alleles as obs_homo2 for pHWE calculation
        for classification in subg_cntl_female[k]:
            if classification == 'obs_homo1':
                phwe_vals[k][0] = sum(subg_cntl_female[k][classification].values())
            elif classification == 'obs_het':
                for key,val in subg_cntl_female[k][classification].items():
                    if not 0 in key[0]:
                        phwe_vals[k][2] += subg_cntl_female[k][classification][key]
                    else:
                        phwe_vals[k][1] += subg_cntl_female[k][classification][key]
            else:
                phwe_vals[k][2] += sum(subg_cntl_female[k][classification].values())

        scores['nClean_' + k] = ",".join((str(x) for x in subg_male[k]['obs_homo1'].values())) +  "," + ",".join((str(x) for x in subg_male[k]['obs_homo2'].values())) +  "," + ",".join((str(x) for x in nclean_female_subg[k])) + ';' + ",".join((str(x) for x in subg_cntl_male[k]['obs_homo1'].values())) + "," + ",".join((str(x) for x in subg_cntl_male[k]['obs_homo2'].values())) + "," + ",".join((str(x) for x in nclean_female_subg_cntl[k]))
        scores['Zhet_' + k] = calc_ExcessHet_multiallelic(zhet_val,zhet_count)
        scores['pHWE_' + k] = calc_pHWE(*phwe_vals[k]) if sum(phwe_vals[k]) >= 5 else '.'
        if type(scores['Zhet_' + k]) == float:
            scores['Zhet_' + k] = "{0:.6f}".format(scores['Zhet_' + k])
        if type(scores['pHWE_' + k]) == float:
            scores['pHWE_' + k] = "{0:.6f}".format(scores['pHWE_' + k])
    return scores


def find_s_d_multiallelic(clean_passing_d, samples): #het_gts):
    
    total_obs_homo_ref = list(clean_passing_d['obs_homo1'].values())[0]
    total_obs_het = sum(list(clean_passing_d['obs_het'].values()))
    total_obs_homo_alt = sum(list(clean_passing_d['obs_homo2'].values()))
    total_obs = total_obs_homo_ref + total_obs_het + total_obs_homo_alt

    maf = 0

    if total_obs > 0:
        maf = (total_obs_het + 2 * total_obs_homo_alt) / (2 * total_obs)
    if maf <= 0.5:
        if total_obs_het == 1 and total_obs_homo_alt == 0:
            idv = find_singleton_multiallelic(samples,clean_passing_d['obs_het'].keys())
            if idv:
                mi.sa.sa_collection[idv].tallySA['singleton'] += 1
        elif total_obs_homo_alt == 1 and total_obs_het == 0:
            idv = find_private_doubleton_multiallelic(samples,clean_passing_d['obs_homo2'].keys())
            if idv: mi.sa.sa_collection[idv].tallySA['p_dblton'] +=1
        elif total_obs_het == 2 and total_obs_homo_alt == 0:
            dbltons = find_doubletons_multiallelic(samples, clean_passing_d['obs_het'].keys())
            for idv in dbltons:
                mi.sa.sa_collection[idv].tallySA['doubleton'] +=1
    else:
        if total_obs_het == 1 and total_obs_homo_ref == 0:
            idv = find_singleton_multiallelic(samples, clean_passing_d['obs_het'].keys())
            if idv:
                mi.sa.sa_collection[idv].tallySA['singleton'] += 1
        elif total_obs_homo_ref == 1 and total_obs_het == 0:
            idv = find_private_doubleton_multiallelic(samples, clean_passing_d['obs_homo1'].keys())
            if idv: mi.sa.sa_collection[idv].tallySA['p_dblton'] += 1
        elif total_obs_het == 2 and total_obs_homo_ref == 0:
            dbltons = find_doubletons_multiallelic(samples,clean_passing_d['obs_het'].keys())
            for idv in dbltons:
                mi.sa.sa_collection[idv].tallySA['doubleton'] += 1


def find_s_d(total_obs, samples):
    """
    """
    maf = 0
    if sum(total_obs) > 0:
        maf = (total_obs[1] + 2 * total_obs[2]) / (2 * sum(total_obs))
    if maf <= 0.5:
        #singleton
        if total_obs[1] == 1 and total_obs[2] == 0:
            idv = find_singleton(samples)
            if idv: mi.sa.sa_collection[idv].tallySA['singleton'] += 1
        #private_doubleton
        elif total_obs[2] == 1 and total_obs[1] == 0:
            idv = find_private_doubleton(samples, 1)
            if idv: mi.sa.sa_collection[idv].tallySA['p_dblton'] += 1
        #doubleton
        elif total_obs[1] == 2 and total_obs[2] == 0:
            dbltons = find_doubletons(samples)
            for idv in dbltons:
                mi.sa.sa_collection[idv].tallySA['doubleton'] += 1
    else:
        #singleton
        if total_obs[1] == 1 and total_obs[0] == 0:
            idv = find_singleton(samples)
            if idv:
                mi.sa.sa_collection[idv].tallySA['singleton'] += 1
        #private_doubleton
        elif total_obs[0] == 1 and total_obs[1] == 0:
            idv = find_private_doubleton(samples, 0)
            if idv: mi.sa.sa_collection[idv].tallySA['p_dblton'] += 1
        #doubleton
        elif total_obs[1] == 2 and total_obs[0] == 0:
            dbltons = find_doubletons(samples)
            for idv in dbltons:
                mi.sa.sa_collection[idv].tallySA['doubleton'] += 1


def find_s_d_chrx(clean_d, samples):
    """
    1/1 Male GT's treated as being 'Heterozygous' since even though males are hemizygous at Chrx, they carry an alt allele
    """
    abc_order = OrderedDict(sorted(mi.sa.sa_collection.items()))      #Used to check SEX of sample 
    maf = 0
    
    total_obs = sum(list(clean_d['male']['obs_homo1'].values())) + sum(list(clean_d['female']['obs_homo1'].values())) + sum(list(clean_d['female']['obs_het'].values())) + sum(list(clean_d['male']['obs_homo2'].values())) + sum(list(clean_d['female']['obs_homo2'].values()))
    total_obs_homo_ref =  sum(list(clean_d['male']['obs_homo1'].values())) + sum(list(clean_d['female']['obs_homo1'].values()))
    total_obs_het = sum(list(clean_d['female']['obs_het'].values())) + sum(list(clean_d['male']['obs_homo2'].values()))
    total_obs_homo_alt =sum(list(clean_d['female']['obs_homo2'].values()))
    female_obs = sum(list(clean_d['female']['obs_homo1'].values())) + sum(list(clean_d['female']['obs_het'].values())) + sum(list(clean_d['female']['obs_homo2'].values()))
    male_obs = sum(list(clean_d['male']['obs_homo1'].values())) + sum(list(clean_d['male']['obs_homo2'].values()))  

    if total_obs > 0:
        maf = (total_obs_het + (2 * total_obs_homo_alt)) / ((2 * female_obs) + male_obs)
    if maf <= 0.5:
        # singleton
        if total_obs_het == 1 and total_obs_homo_alt == 0:
            if sum(list(clean_d['female']['obs_het'].values())):                    #Het sample is female, use female GT 0/1
                idv = find_singleton_chrx(samples, clean_d['female']['obs_het'].keys(), "1")
                if idv: mi.sa.sa_collection[idv].tallySA['singleton'] += 1
            if sum(list(clean_d['male']['obs_homo2'].values())):                  #Het sample is male, use male GT 1/1
                idv = find_singleton_chrx(samples, clean_d['male']['obs_homo2'].keys(), "0")
                if idv: mi.sa.sa_collection[idv].tallySA['singleton'] += 1
         #private_doubleton
        elif total_obs_homo_alt == 1 and total_obs_het == 0:
            idv = find_private_doubleton_chrx(samples,clean_d['female']['obs_homo2'].keys(), "1")
            if idv: mi.sa.sa_collection[idv].tallySA['p_dblton'] +=1

        #doubleton
        elif total_obs_het == 2 and total_obs_homo_alt == 0:
            if sum(list(clean_d['male']['obs_homo2'].values())):                #If Male hets, use 1/1 GT
                dbltons = find_doubletons_chrx(samples,clean_d['male']['obs_homo2'].keys(), "0")
                for idv in dbltons:
                    mi.sa.sa_collection[idv].tallySA['doubleton'] += 1
            if sum(list(clean_d['female']['obs_het'].values())):                #If Female hets, use 0/1 GT
                dbltons = find_doubletons_chrx(samples,clean_d['female']['obs_het'].keys(), "1")
                for idv in dbltons:
                    mi.sa.sa_collection[idv].tallySA['doubleton'] += 1

    else:
        #singleton
        if total_obs_het == 1 and total_obs_homo_ref == 0:
            if sum(list(clean_d['female']['obs_het'].values())):                #If Female hets, use 0/1 GT
                idv = find_singleton_chrx(samples, clean_d['female']['obs_het'].keys(), "1")
                if idv: mi.sa.sa_collection[idv].tallySA['singleton'] += 1
            if sum(list(clean_d['male']['obs_homo2'].values())):              #If Male hets, use 1/1 GT
                idv = find_singleton_chrx(samples, clean_d['male']['obs_homo2'].keys(), "0")
                if idv: mi.sa.sa_collection[idv].tallySA['singleton'] += 1

        #private_doubleton
        elif total_obs_homo_ref == 1 and total_obs_het == 0:
            if sum(list(clean_d['male']['obs_homo1'].values())):                #If male Homozygous_ref, use 0/0 GT
                idv = find_private_doubleton_chrx(samples, clean_d['male']['obs_homo1'].keys(), "0")
                if idv: mi.sa.sa_collection[idv].tallySA['p_dblton'] += 1
            else:                                                               #If female Homozygous_ref, use 0/0 GT
                idv = find_private_doubleton_chrx(samples, clean_d['female']['obs_homo1'].keys(), "1")
                if idv: mi.sa.sa_collection[idv].tallySA['p_dblton'] += 1

        #doubleton
        elif total_obs_het == 2 and total_obs_homo_ref == 0:                    #If Male Het, use 1/1 GT
            if sum(list(clean_d['male']['obs_homo2'].values())):
                dbltons = find_doubletons_chrx(samples,clean_d['male']['obs_homo2'].keys(), "0")
                for idv in dbltons:
                    mi.sa.sa_collection[idv].tallySA['doubleton'] += 1
            if sum(list(clean_d['female']['obs_het'].values())):                #If Female het, use 0/1 GT
                dbltons = find_doubletons_chrx(samples,clean_d['female']['obs_het'].keys(), "1")
                for idv in dbltons:
                    mi.sa.sa_collection[idv].tallySA['doubleton'] += 1

def find_singleton(samples):
    for k, sm in samples.items():
        if sm['GT'] in {(0,1), (1,0)}:
            try:
                if (sm['DP'] >= cfg.MINDP
                   and sm['GQ'] >= cfg.MINGQ):
                        return k
            except TypeError:  # TypeError: unorderable types: NoneType() < int() (missing DP)
                continue

def find_singleton_multiallelic(samples, het_alleles):
    for k, sm in samples.items():
        for gts in het_alleles:
            if sm['GT'] in gts:
                try:
                    if (sm['DP'] >= cfg.MINDP
                        and sm['GQ'] >= cfg.MINGQ):
                            return k
                except TypeError:  # TypeError: unorderable types: NoneType() < int() (missing DP)
                    continue

def find_singleton_chrx(samples, het_alleles, gender):
    abc_order = OrderedDict(sorted(mi.sa.sa_collection.items()))
    for k, sm in samples.items():
        if gender == abc_order[k].details_dict.SEX:
            for gts in het_alleles:
                if sm['GT'] in gts:
                    try:
                        if (sm['DP'] >= cfg.MINDP
                            and sm['GQ'] >= cfg.MINGQ):
                                return k
                    except TypeError:  # TypeError: unorderable types: NoneType() < int() (missing DP)
                        continue

def find_private_doubleton(samples, allele):
    for k, sm in samples.items():
        if sm['GT'] == (allele, allele):
            try:
                if (sm['DP'] >= cfg.MINDP
                   and sm['GQ'] >= cfg.MINGQ):
                    return k
            except TypeError:  # TypeError: unorderable types: NoneType() < int() (missing DP)
                continue

def find_private_doubleton_multiallelic(samples, allele):
    for k, sm in samples.items():
        for gts in allele:
            if (sm['GT']) in gts:
                try:
                    if (sm['DP'] >= cfg.MINDP
                        and sm['GQ'] >= cfg.MINGQ):
                        return k
                except TypeError:  # TypeError: unorderable types: NoneType() < int() (missing DP)
                    continue


def find_private_doubleton_chrx(samples, allele, gender):
    abc_order = OrderedDict(sorted(mi.sa.sa_collection.items()))
    for k, sm in samples.items():
        for gts in allele:
            if (sm['GT']) in gts:
                try:
                    if gender == abc_order[k].details_dict.SEX:
                        if (sm['DP'] >= cfg.MINDP
                            and sm['GQ'] >= cfg.MINGQ):
                            return k
                except TypeError:  # TypeError: unorderable types: NoneType() < int() (missing DP)
                    continue

def find_singleton_chrx(samples, het_alleles, gender):
    abc_order = OrderedDict(sorted(mi.sa.sa_collection.items()))
    for k, sm in samples.items():
        if gender == abc_order[k].details_dict.SEX:
            for gts in het_alleles:
                if sm['GT'] in gts:
                    try:
                        if (sm['DP'] >= cfg.MINDP
                            and sm['GQ'] >= cfg.MINGQ):
                                if gender == abc_order[k].details_dict.SEX:
                                    return k
                    except TypeError:  # TypeError: unorderable types: NoneType() < int() (missing DP)
                        continue

def find_doubletons(samples):
    k_list = list()
    for k,sm in samples.items():
        if sm['GT'] in {(0,1), (1,0)}:
            try:
                if (sm['DP'] >= cfg.MINDP
                   and sm['GQ'] >= cfg.MINGQ):
                    k_list.append(k)
            except TypeError:  #TypeError: unorderable types: NoneType() < int() (missing DP)
                continue

        if len(k_list) == 2:
            return k_list

    return k_list
 
def find_doubletons_multiallelic(samples, het_alleles):
    k_list = list()
    for k,sm in samples.items():
        for gts in het_alleles:
            if sm['GT'] in gts:
                try:
                    if (sm['DP'] >= cfg.MINDP
                        and sm['GQ'] >= cfg.MINGQ):
                            k_list.append(k)
                except TypeError:  #TypeError: unorderable types: NoneType() < int() (missing DP)
                    continue
            
            if len(k_list) == 2:
                return k_list

    return k_list

def find_doubletons_chrx(samples, het_alleles, gender):
    abc_order = OrderedDict(sorted(mi.sa.sa_collection.items()))
    k_list = list()
    for k,sm in samples.items():
        if gender == abc_order[k].details_dict.SEX:
            for gts in het_alleles:
                if sm['GT'] in gts:
                    try:
                        if (sm['DP'] >= cfg.MINDP
                            and sm['GQ'] >= cfg.MINGQ):
                                k_list.append(k)
                    except TypeError:  #TypeError: unorderable types: NoneType() < int() (missing DP)
                        continue
                if len(k_list) == 2:
                    return k_list

    return k_list

def gather_intersect_fam_vcf_samples(vcf_samples, fam_samples):
    """
    """
    # Method 5 - set, one pass
    for key,sm in vcf_samples.items():
        for subset, sm_set in fam_samples.items():
            if key in sm_set['set']:
                fam_samples[subset]['dict'][key] = sm
    return fam_samples


def check_mendelian_errors_multiallelic(prefix, rec):
    """
    If 1 parent:
        If child_allele 1 or child_allele 2 != parent_allele 1 or parent_allele 2: Mendelian error

    If 2 parents:
        If child_allele 1 or child_allele 2 != parent1_allele1 or parent1_allele2: Mendelian error
        If child_allele 1 or child_allele 2 != parent2_allele1 or parent2_allele2: Mendelian error

    """
    samples = rec.samples
    mend_pairs = 0
    mend_error = 0
    for kid in mi.sa.get_mi_kids():
        validparents = []
        found_mi_error = 0
        father = mi.sa.get_father(kid)
        mother = mi.sa.get_mother(kid)
        if father in mi.sa.id_list and is_good_gt(samples[ father ]):
            mend_pairs += 1
            validparents.append(father)
        if mother in mi.sa.id_list and is_good_gt(samples[ mother ]):
            mend_pairs += 1
            validparents.append(mother)
            #print('valid parent found')
        # mend_error is child having allele not from parents
        if validparents == []:
            continue
        elif len(validparents) == 1: # Only 1 parent
            kid_allele1 = samples[kid]['GT'][0]
            kid_allele2 = samples[kid]['GT'][1]

            if ( kid_allele1 == samples[validparents[0]]['GT'][0] ) or ( kid_allele1 == samples[validparents[0]]['GT'][1] ): #Match
                continue
            elif ( kid_allele2 == samples[validparents[0]]['GT'][0] ) or ( kid_allele2 == samples[validparents[0]]['GT'][1] ):
                continue
            else:
                mend_error += 1
                found_mi_error = 1
        elif len(validparents) == 2: # 2 parents
            kid_allele1 = samples[kid]['GT'][0]
            kid_allele2 = samples[kid]['GT'][1]

            if ( kid_allele1 == samples[validparents[0]]['GT'][0] ) or ( kid_allele1 == samples[validparents[0]]['GT'][1] ): #Kid allele 1 matches parent 1 allele
                if ( kid_allele2 == samples[validparents[1]]['GT'][0] ) or ( kid_allele2 == samples[validparents[1]]['GT'][1] ): #Kid allele 2 matches parent 2 allele
                    continue # Clean - Both kid alleles match a parent allele
                else:
                    pass #Kid allele 1 matches  parent_1 allele, but kid allele 2 doesn't match parent_2 allele. 
            if ( kid_allele2 == samples[validparents[0]]['GT'][0] ) or ( kid_allele2 == samples[validparents[0]]['GT'][1]) : # Kid allele 2 matches parent 1 allele
                if ( kid_allele1 == samples[validparents[1]]['GT'][0] ) or ( kid_allele1 == samples[validparents[1]]['GT'][1] ): #Kid allele 1 matches parent 2 allele
                    continue # Clean - Both kid alleles match a parent allele
                else: # Mendelian Error - only 1 kid allele matches a parent allele
                    mend_error += 1
                    found_mi_error = 1
            else: # Mendelian Error - no kid alleles match either parent alleles
                mend_error += 1
                found_mi_error = 1

            #else: clean child

        mi.sa.sa_collection[kid].tallySA['mend_pair'] += 1
        if found_mi_error:
            # save details
            p2 = samples[ validparents[1] ]['GT'] if len(validparents)>1 else ('.')
            write_mendelian_errors(prefix, rec,
                                   [mi.sa.get_fam_id(kid), validparents[0], kid],
                                   [samples[ validparents[0] ]['GT'], p2, samples[kid]['GT'] ]
                                   )
            mi.sa.sa_collection[kid].tallySA['vp' + str(len(validparents)) ] += 1
    return mend_pairs, mend_error
#BI ALLELIC


def check_mendelian_errors(prefix, rec):
    """
    """
    abc_order = OrderedDict(sorted(mi.sa.sa_collection.items()))
    samples = rec.samples
    mend_pairs = 0
    mend_error = 0
    for kid in mi.sa.get_mi_kids():
        validparents = []
        found_mi_error = 0

        father = mi.sa.get_father(kid)
        mother = mi.sa.get_mother(kid)


        if father in mi.sa.id_list and is_good_gt(samples[ father ]):
            mend_pairs += 1
            validparents.append(father)

        if mother in mi.sa.id_list and is_good_gt(samples[ mother ]):
            mend_pairs += 1
            validparents.append(mother)
        
        # mend_error is child having allele not from parents
        if validparents == []:
            continue
        elif len(validparents) == 1:
            if (samples[kid]['GT'] in {(0,0), (1,1)}
                and
                samples[ validparents[0] ]['GT'] in {(0,0), (1,1)}
                and
                (abs(sum(samples[kid]['GT']) - sum(samples[ validparents[0] ]['GT'])) == 2)
                ):
                mend_error += 1
                found_mi_error = 1
        elif len(validparents) == 2:
            # Two parent mendelian errors
            # Case 1) 0/0 & 1/1: no 1/1, 0/0
            # Case 2) 0/0 & 0/0: no 0/1, 1/1
            # Case 3) 0/0 & 0/1: no 1/1
            # Case 4) 1/1 & 1/1: no 0/1, 0/0
            # Case 5) 0/1 & 1/1: no 0/0
            sum_p_genos = sum(samples[ validparents[0] ]['GT']) + sum(samples[ validparents[1] ]['GT'])
            sum_k_geno = sum(samples[kid]['GT'])

            # case #2, 3
            if sum_p_genos in [0, 1] and sum_k_geno > sum_p_genos:
                mend_error += 1
                found_mi_error = 1
            # case #4
            elif sum_p_genos == 4 and sum_k_geno in [0, 1]:
                mend_error += 1
                found_mi_error = 1
            # case #5
            elif sum_p_genos == 3 and sum_k_geno == 0:
                mend_error += 1
                found_mi_error = 1
            # case #1
            elif (samples[kid]['GT'] in {(0,0), (1,1)}
                and samples[ validparents[0] ]['GT'] != samples[ validparents[1] ]['GT']
                and sum_p_genos == 2):
                mend_error += 1
                found_mi_error = 1

            #else: clean child
        
        mi.sa.sa_collection[kid].tallySA['mend_pair'] += 1

        if found_mi_error:
            p2 = samples[ validparents[1] ]['GT'] if len(validparents)>1 else ('.')

            write_mendelian_errors(prefix, rec,
                                   [mi.sa.get_fam_id(kid), validparents[0], kid],
                                   [samples[ validparents[0] ]['GT'], p2, samples[kid]['GT'] ]
                                   )
            mi.sa.sa_collection[kid].tallySA['vp' + str(len(validparents)) ] += 1
    return mend_pairs, mend_error

def check_mendelian_errors_chrx(prefix, rec):
    """
    If 1 parent:
        If child_allele 1 or child_allele 2 != parent_allele 1 or parent_allele 2: Mendelian error

    If 2 parents:
        If child_allele 1 or child_allele 2 != parent1_allele1 or parent2_allele2: Mendelian error
        If child_allele 1 or child_allele 2 != parent2_allele1 or parent2_allele1: Mendelian error

    """
    abc_order = OrderedDict(sorted(mi.sa.sa_collection.items()))
    samples = rec.samples
    mend_pairs = 0
    mend_error = 0
    for kid in mi.sa.get_mi_kids():
        if abc_order[kid].details_dict.SEX == "0":
            if not samples[kid]['GT'][0] == samples[kid]['GT'][1]:
                continue
        validparents = []
        found_mi_error = 0
        father = mi.sa.get_father(kid)
        mother = mi.sa.get_mother(kid)
        if father in mi.sa.id_list and is_good_gt(samples[ father ]):
            father_allele1 = samples[father]['GT'][0]
            father_allele2 = samples[father]['GT'][1]
            if father_allele1 == father_allele2:
                if abc_order[kid].details_dict.SEX == "1": #It is female
                    mend_pairs += 1
                    validparents.append(father)
        if mother in mi.sa.id_list and is_good_gt(samples[ mother ]):
            mend_pairs += 1
            validparents.append(mother)
        if validparents == []:
            continue
        elif len(validparents) == 1: # Only 1 parent
            kid_allele1 = samples[kid]['GT'][0]
            kid_allele2 = samples[kid]['GT'][1]
            if ( kid_allele1 == samples[validparents[0]]['GT'][0] ) or ( kid_allele1 == samples[validparents[0]]['GT'][1] ): #Match
                continue
            elif ( kid_allele2 == samples[validparents[0]]['GT'][0] ) or ( kid_allele2 == samples[validparents[0]]['GT'][1] ):
                continue
            else:
                #print(f' MENDELIAN ERROR')
                mend_error += 1
                found_mi_error = 1
        elif len(validparents) == 2: # 2 parents
            kid_allele1 = samples[kid]['GT'][0]
            kid_allele2 = samples[kid]['GT'][1]
            if ( kid_allele1 == samples[validparents[0]]['GT'][0] ) or ( kid_allele1 == samples[validparents[0]]['GT'][1] ): #Kid allele 1 matches parent 1 allele
                if ( kid_allele2 == samples[validparents[1]]['GT'][0] ) or ( kid_allele2 == samples[validparents[1]]['GT'][1] ): #Kid allele 2 matches parent 2 allele
                    continue # Clean - Both kid alleles match a parent allele
                else:
                    pass #Kid allele 1 matches  parent_1 allele, but kid allele 2 doesn't match parent_2 allele. 
            if ( kid_allele2 == samples[validparents[0]]['GT'][0] ) or ( kid_allele2 == samples[validparents[0]]['GT'][1]) : # Kid allele 2 matches parent 1 allele
                if ( kid_allele1 == samples[validparents[1]]['GT'][0] ) or ( kid_allele1 == samples[validparents[1]]['GT'][1] ): #Kid allele 1 matches parent 2 allele
                    continue # Clean - Both kid alleles match a parent allele
                else: # Mendelian Error - only 1 kid allele matches a parent allele
                   # print(f' MENDELIAN ERROR 2')
                    mend_error += 1
                    found_mi_error = 1
            else: # Mendelian Error - no kid alleles match either parent alleles
                mend_error += 1
                found_mi_error = 1

            #else: clean child

        mi.sa.sa_collection[kid].tallySA['mend_pair'] += 1
        if found_mi_error:
            
            # save details
            p2 = samples[ validparents[1] ]['GT'] if len(validparents)>1 else ('.')
            write_mendelian_errors(prefix, rec,
                                   [mi.sa.get_fam_id(kid), validparents[0], kid],
                                   [samples[ validparents[0] ]['GT'], p2, samples[kid]['GT'] ]
                                   )
            mi.sa.sa_collection[kid].tallySA['vp' + str(len(validparents)) ] += 1
    return mend_pairs, mend_error

if __name__ == "__main__":
    main()
    # cProfile.run('main()', None, 'tottime')
