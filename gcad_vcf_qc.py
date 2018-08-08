#!/usr/bin/env python3
#-*- coding: utf-8 -*-

from argparse import ArgumentParser
from pysam import VariantFile

#import configparser
import csv
import os
from collections import namedtuple, OrderedDict, Counter
import time
import cProfile

from utils.VariantAnnotation.vflags import calcVA
from utils.stats.count_gt import is_good_gt

import utils.SampleAnnotation.sample_annotation as mi
import config as cfg

import warnings
warnings.simplefilter('always')



def extractSubSets(fam):
    """
    extractSubSets - Creates a OrderedDict() where Subset groupings are the key, 
                     values are SampID within it
    @return Samples Dict per Subset
    """
    SampleFamDetail = namedtuple('SampleFamDetail',['SampID','FID','SubjID','FA','MO',
                                                    'SEX','AFF','AD','AGE','ADSPWGS',
                                                    'Subset','Subgroup','Race_Ethnicity','SeqCtr','ExcludeFromZHet'])

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

def write_subset_stats(subset, rec, vf, abhet, passing, failing, missing, gt_failed, depth_sum, clean_obs, mend_pairs, mend_errors):
    """
    """
    outfile = 'summary.snv.{}.tsv'.format( subset)
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
        mean_depth = depth_sum / sum_clean if sum_clean else 0

        writer.writerow({'CHR': rec.contig, 'POS': rec.pos,
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
                         'QUAL':"{0:.2f}".format(rec.qual),'FILTER':",".join(rec.filter.keys()),
                         })


    return

def write_mendelian_errors(rec, fam_info, genos ): # mmmm, genos
    """
    """
    outfile = 'summary.mi.tsv'
    newfile = not os.path.exists(outfile)

    with open(outfile, 'a') as csvfile:
        fieldnames = ['FID','PID','CID','CHR','POS','REF','ALT','PGT','CGT',]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames , delimiter='\t', lineterminator='\n')

        if newfile:
            writer.writeheader()

        writer.writerow({'FID': fam_info[0], 'PID':fam_info[1], 'CID': fam_info[2],
                         'CHR': rec.contig, 'POS': rec.pos,
                         'REF': rec.ref, 'ALT':rec.alts[0],
                         'PGT': "/".join(map(str,genos[0])),'CGT': "/".join(map(str,genos[1]))
                         })

def write_indiv_summary():
    """
    """
    outfile = 'summary.indiv.tsv'

    with open(outfile, 'w') as csvfile:
        fieldnames = ['SampleID','SEX',
                      'total_nRR','total_nRA','total_nAA','Missing','Set_Missing',
                      'Singleton','Doubleton','HetHom',
                      'Ti','Tv','TiTvRatio','IndMeanDepth',
                      '1P_MI','2P_MI','MI_pairs',]
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
            mean_depth = val.dp_total / good_gt if good_gt else 0

            writer.writerow({'SampleID': indiv, 'SEX':val.details_dict.SEX,
                            'total_nRR': val.tallySA[(0,0)],'total_nRA': good_het_gt,'total_nAA': val.tallySA[(1,1)],
                            'Missing': val.tallySA[(None,None)],'Set_Missing': val.tallySA[-9],
                            'Singleton': val.tallySA['singleton'],
                            #'Private_Doubleton': val.tallySA['p_dblton'],
                            'Doubleton': val.tallySA['doubleton'] + val.tallySA['p_dblton'],
                            'HetHom':"{0:.2f}".format(het_hom),
                            'Ti':val.tallySA['ti'], 'Tv':val.tallySA['tv'], 'TiTvRatio':"{0:.2f}".format(ti_tv),'IndMeanDepth':"{0:.2f}".format(mean_depth),
                            '1P_MI':val.tallySA['vp1'],'2P_MI':val.tallySA['vp2'],'MI_pairs':val.tallySA['mend_pair'],
                            })


def delete_previous_outputs(prefix, subsets):
    """
    """
    mi_file = 'summary.mi.tsv'
    for out_file in [mi_file] + [ '{}.{}.tsv'.format(prefix, x) for x in subsets]:
        if os.path.exists(out_file):
            os.remove(out_file)
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
    grp_settings.add_argument('--hetz_lim1', type=int, help='', default=99999, required=False)
    grp_settings.add_argument('--hetz_lim2', type=int, help='', default=99999, required=False)
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

    mi.createSampleAnnotation(args.fam)
    samplesDict, famCt = extractSubSets(args.fam) # returns dict

    print("[FAM] Found {} subsets: {}; totaling {} sampIDs".format(len(samplesDict.keys()), list(samplesDict.keys()), famCt))
    print("[FAM] {}".format([  "{}:{}".format(k, len(samplesDict[k]))  for k in samplesDict.keys()]))
    #{k:rec.samples[k] for k in samples['sub1']}
    #{k:rec.samples[k] for k in samples['sub1'] if k in rec.samples}
    #{key: d[key] for key in d.viewkeys() & l}
    print("[FAM] {} unique subgroups:{}".format(len(mi.sa.subgroups),sorted( mi.sa.subgroups)))
    print("[FAM] {} ".format([  "{}:{}".format(k, v)  for k,v in mi.sa.subsets.items() ]))

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

    delete_previous_outputs('summary.snv',list(samplesDict.keys()))

    for rec in vcf_in.fetch(rChr, rStart, rEnd):
        #vcf_out.write(rec)
        #if ct>500:break
        #if rec.pos < 10684424: continue
        if len(rec.alts) > 1:
            print("Warning found multiallelic variant")
            contine

        #print(str(rec.contig) + '\t', str(rec.pos) + '\t', str(rec.ref) + '\t', str(rec.alts[0]) + '\t', end='')

        #
        samplesDict = gather_intersect_fam_vcf_samples(rec.samples, samplesDict)
        grp_obs = list()

        for subset, sm_list in samplesDict.items():
            [vf, abhet, passing, failing, missing, gt_failed, depth_sum, clean_obs, subg] = calcVA(sm_list['dict'], rec)

            grp_obs.append(clean_obs)
            mend_pairs, mend_errors = check_mendelian_errors(rec)
            #subg = calculate_subgroup_scores(subg)

            #if subset == 'ADSPfamWGS' and mend_errors > 0:
            #    print("VFLAGS_{}={};".format(subset, vf), end='')
            #    print("ABHet_{}={};".format(subset, abhet), end='')
                #print(str(rec.contig) + '\t', str(rec.pos) + '\t', str(rec.ref) + '\t', str(rec.alts[0]) + '\t', end='')
                #print("failing={};gt_failed={}".format(failing, gt_failed), end=' ')
                #print("mi={};mp={}".format(mend_errors, mend_pairs), end=' ')
                #print()

            write_subset_stats(subset, rec, vf, abhet, passing, failing, missing, gt_failed, depth_sum, clean_obs, mend_pairs, mend_errors )

        #if len(mi.sa.singletons) > 1 :
            #for idv in mi.sa.singletons:
                #mi.sa.sa_collection[idv].tallySA['singleton'] -= 1
            #mi.sa.singletons.clear()
        #if len(mi.sa.private_dbltons) > 1:
            #for idv in mi.sa.private_dbltons:
                #mi.sa.sa_collection[idv].tallySA['p_dblton'] -= 1
            #mi.sa.private_dbltons.clear()
        #if len(mi.sa.dbltons) > 2:
            #for idv in mi.sa.dbltons:
                #mi.sa.sa_collection[idv].tallySA['doubleton'] -= 1
            #mi.sa.dbltons.clear()

        total_obs = list(map(sum, zip(*grp_obs)))

        find_s_d(total_obs, rec.samples)

        #print()
        ct += 1

    end = time.time()
    print("{0:.2f}".format(end - start))
    write_indiv_summary()

def find_s_d(total_obs, samples):
    maf = 0
    if sum(total_obs) > 0:
        maf = (total_obs[1] + 2 * total_obs[2]) / (2 * sum(total_obs))
    if maf <= 0.5:
        #singleton
        if total_obs[1] == 1 and total_obs[2] == 0:
            idv = find_singleton(samples)
            if idv: mi.sa.sa_collection[idv].tallySA['singleton'] += 1
        elif total_obs[2] == 1 and total_obs[1] == 0:
            idv = find_private_doubleton(samples, 1)
            if idv: mi.sa.sa_collection[idv].tallySA['p_dblton'] += 1
        elif total_obs[1] == 2 and total_obs[2] == 0:
            dbltons = find_doubletons(samples)
            for idv in dbltons:
                mi.sa.sa_collection[idv].tallySA['doubleton'] += 1
    else:
        if total_obs[1] == 1 and total_obs[0] == 0:
            idv = find_singleton(samples)
            mi.sa.sa_collection[idv].tallySA['singleton'] += 1
        elif total_obs[0] == 1 and total_obs[1] == 0:
            idv = find_private_doubleton(samples, 0)
            if idv: mi.sa.sa_collection[idv].tallySA['p_dblton'] += 1
        elif total_obs[1] == 2 and total_obs[0] == 0:
            dbltons = find_doubletons(samples)
            for idv in dbltons:
                mi.sa.sa_collection[idv].tallySA['doubleton'] += 1

def find_singleton(samples):
    for k,sm in samples.items():
        if (sm['GT'] == (1, 0)
            or sm['GT'] == (0, 1)
            and (sm['DP'] >= cfg.MINDP
                and sm['GQ'] >= cfg.MINGQ)
            ):
            return k

def find_private_doubleton(samples, allele):
    for k,sm in samples.items():
        if (sm['GT'] == (allele, allele)
            and (sm['DP'] >= cfg.MINDP
                and sm['GQ'] >= cfg.MINGQ)):
            return k

def find_doubletons(samples):
    k_list = list()
    for k,sm in samples.items():
        if (sm['GT'] == (1, 0)
            or sm['GT'] == (0, 1)
            and (sm['DP'] >= cfg.MINDP
                and sm['GQ'] >= cfg.MINGQ)
            ):
            k_list.append(k)
        if len(k_list) == 2:
            return k_list

    return k_list
def gather_intersect_fam_vcf_samples(vcf_samples, fam_samples):
    """
    """
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
    for key,sm in vcf_samples.items():
        for subset, sm_set in fam_samples.items():
            if key in sm_set['set']:
                fam_samples[subset]['dict'][key] = sm
    return fam_samples

def check_mendelian_errors(rec):
    """
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
            # mend error when homozygous child that doesn't allele match homozygous parents
            if (samples[kid]['GT'] in {(0,0), (1,1)}
                and
                (sum(samples[ validparents[0] ]['GT']) + sum(samples[ validparents[1] ]['GT']) in {0,4} )
                and
                (samples[kid]['GT'] != samples[ validparents[0] ]['GT'])
                ):
                mend_error += 1
                found_mi_error = 1
            elif sum(samples[kid]['GT']) == 1:
            # mend_error if heterozygous child has homozygous matching parents
                if sum(samples[father]['GT']) + sum(samples[mother]['GT']) in {0,4}:
                    mend_error += 1
                    found_mi_error = 1
            #else: clean child

        mi.sa.sa_collection[kid].tallySA['mend_pair'] += 1

        if found_mi_error:
            # save details
            write_mendelian_errors(rec,
                                   [mi.sa.get_fam_id(kid), validparents[0], kid],
                                   [samples[ validparents[0] ]['GT'], samples[kid]['GT'] ]
                                   )
            mi.sa.sa_collection[kid].tallySA['vp' + str(len(validparents)) ] += 1

    return mend_pairs, mend_error

if __name__ == "__main__":
    main()
    cProfile.run('main()', None, 'cumtime')
