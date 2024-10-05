#!/usr/bin/env python3

import config as cfg
import utils.SampleAnnotation.sample_annotation as mi
from collections import OrderedDict
from itertools import product
import copy

def count_gt(samples, rec_details, in_region):
    """
    count_gt - sample genotype {0/0, 0/1, 1/1}
               Apply Genotype-level QC: DP<10, GQ<20; set to ./.

    Parameters:
    samples: samples object
    rec_details (dict): record details
    in_region (bool): is the record inside of a WES target region

    Returns:
    obs_hom1 (int): count of observed 0/0
    obs_hts  (int): count of observerd 0/1
    obs_hom2 (int): count of observed 1/1
    missing  (int): count of sample missing a genotype
    gt_failed (int): count of sample's genotype not passing DP/GQ minimum thresholds
    depth_sum (int): sum of read depth for variant, from all passing samples
    failed   (list): count of samples's genotype not passing miniumum DP/GQ, by genotype [0/0, 0/1, 1/1]
    het_ad   (int): allele depth sum for all passing samples
    het_dp   (int): sum depth (DP) for all passing samples
    subgroup_counts (list): genotype count by subgroup
    subgroup_counts_cntrls (list): genotype count by subgroup only on controls

    Also alters global objects: mi.sa.tally

    """
    obs_hts = 0
    obs_hom1 = 0
    obs_hom2 = 0
    missing = 0
    depth_sum = 0
    gt_failed = 0
    failed = [0,0,0]
    het_ad = 0 # for ABHET
    het_dp = 0 # for ABHET

    #subgroup_counts = OrderedDict({key:[0,0,0] for key in mi.sa.subgroups})
    #subgroup_counts_cntrls = copy.deepcopy(subgroup_counts)

    ref = rec_details['ref']
    alt = rec_details['alt'][0]

    mi.sa.clear_mpairs()

    for k,sm in samples.items():


        if (None in sm['GT']
            ):
            missing += 1
            sm['GT'] = (None, None)
            tallyMissingSample(k, sm)
            continue

        try:
            if (
                #sm['DP'] == None
                #or
                sm['DP'] < cfg.MINDP
                or sm['GQ'] < cfg.MINGQ
                ):
                if sm['GT'] == (0,0):
                    failed[0] += 1
                elif (sm['GT'] in {(0,1), (1,0)} ):
                    failed[1] += 1
                elif sm['GT'] == (1,1):
                    failed[2] += 1
                else:
                    raise Exception('Unknown GT in sample: {} {}'.format(sm['GT'], k))

                depth_sum += sm['DP']
                sm['GT'] = (None, None)
                tallyFailedSample(k, sm)
                gt_failed += 1
                continue
        except TypeError:  #TypeError: unorderable types: NoneType() < int() (missing DP)
            #sm['DP'] == None
            if sm['GT'] == (0,0):
                failed[0] += 1
            elif (sm['GT'] in {(0,1), (1,0)} ):
                failed[1] += 1
            elif sm['GT'] == (1,1):
                failed[2] += 1
            else:
                raise Exception('Unknown GT in sample: {} {}'.format(sm['GT'], k))

            #sm['GT'] += (0)
            sm['GT'] = (None, None)
            tallyFailedSample(k, sm)
            gt_failed += 1
            #raise TypeError("Weird {},{}, {}, {}".format(k,str(sm['GT']), str(sm['DP']), str(sm['AD']) ))
            continue
        except:
            raise

        depth_sum += sm['DP']
        tallyIndelSample(k, ref, alt)

        if (sm['GT'] in {(0,1), (1,0)} ):
            obs_hts += 1
            het_ad += sm['AD'][0]
            het_dp += sm['DP']

            #subgroup_counts, subgroup_counts_cntrls = increment_subgroup(k, 1, subgroup_counts, subgroup_counts_cntrls)
            tallyTiTv(k, sm, ref, alt, in_region)
        elif sm['GT'] == (0,0):
            obs_hom1 += 1
            #subgroup_counts, subgroup_counts_cntrls = increment_subgroup(k, 0, subgroup_counts, subgroup_counts_cntrls)
        elif sm['GT'] == (1,1):
            obs_hom2 += 1
            #subgroup_counts, subgroup_counts_cntrls = increment_subgroup(k, 2, subgroup_counts, subgroup_counts_cntrls)
            #if obs_hom2 == 1 and obs_hts == 0:
                #p_dblton = k
            tallyTiTv(k, sm, ref, alt, in_region)
        else:
            raise TypeError("Weird GT")

        tallyPassingSample(k, sm)


    return [obs_hom1, obs_hts, obs_hom2, missing, gt_failed, depth_sum, failed, het_ad, het_dp]

def count_gt_multiallelic(samples,rec_details):
    """
    count_gt_multiallelic - sample genotype {0/0, 0/1, 0/2, 0/3, 0/4, 0/5, 0/6, 1/1, 1/2 ...}
            Apply Genotype-level QC: DP<10, GQ<20; set to ./.

    Parameters:
    samples: samples object
    rec_details (dict): record details

    Returns:
    passing_d      (dict): counts of observed passing Homozygous_reference, Heterozygous, Homozygous_alt Genotypes
    failing_d      (dict): counts of samples's Homozygous_reference, Heterozygous, Homozygous_alt Genotypes not passing minimum DP/GQ 
    missing         (int): count of samples missing a genotype
    gt_failed       (int): count of sample's genotype not passing DP/GQ minimum thresholds
    clean_passing_d (dict): counts of observed passing Homozygous_reference, Heterozygous, Homozygous_alt GTs (Will not be set to 0 in summary.snv output file if vlfags present)
    depth_sum        (int): sum of read depth for variant, from all passing samples
    abhet_AD_list   (list): List of allele depth sum for all passing samples per allele
    abhet_DP_lit    (list): List of sum depth (DP) for all passing samples for all alleles
    subgroup_counts (list): genotype count by subgroup
    subgroup_counts_cntrls (list): genotype count by subgroup only on controls
    allele_count_dict          (dict): counts of alleles of subjects and controls
    zhet_sample_counts (dict): Het and Homozygous_alt counts of subjects and controls

    Also alters global objects: mi.sa.tally
    """
    N = len(rec_details['alt'])
    allele_list = [n for n in range(0,N+1)]
    missing = 0
    depth_sum = 0
    gt_failed = 0
    #zhet_sample_counts = OrderedDict({key:[0,0] for key in mi.sa.subgroups}) # Counts for Het and Homozygous_alts
    abhet_AD_list = [0 for i in allele_list]
    abhet_DP_list = copy.deepcopy(abhet_AD_list)
    #subgroup_counts = OrderedDict({key:{'obs_homo1':{},'obs_het':{},'obs_homo2':{}} for key in mi.sa.subgroups})
    #subgroup_counts_cntrls = copy.deepcopy(subgroup_counts)
    mi.sa.clear_mpairs()
    clean_passing_d = {'obs_homo1':{},'obs_het':{},'obs_homo2':{}}
    #allele_count_dict = OrderedDict({key:[0 for n in range(0,N+1)] for key in mi.sa.subgroups})

    passing_d = {'obs_homo1': {((i, b), (b, i)): 0 for i in allele_list for b in allele_list[i:] if (i, b) == (0, 0) and (b, i) == (0, 0)},
                            'obs_het': {((i, b), (b, i)): 0 for i in allele_list for b in allele_list[i:] if (i, b) != (b, i)},
                            'obs_homo2': {((i, b), (b, i)): 0 for i in allele_list for b in allele_list[i:] if (i, b) == (b, i) and i != 0}}
    failing_d = copy.deepcopy(passing_d)
    
    for i in allele_list:    #Create Passing and Failing Dictionary with keys being all possible Genotypes. Pipeline coded to recognize both orientations of GT's example: (0,1) and (1,0) GT's
            for b in allele_list[i:]:
                for subgroup in mi.sa.subgroups:
                    gt_category = ('obs_homo1' if (i, b) == (0, 0) and (b, i) == (0, 0) 
                    else 'obs_het' if (i, b) != (b, i) 
                    else 'obs_homo2')
                    
                    #subgroup_counts[subgroup][gt_category][(i, b), (b, i)] = 0
                    #subgroup_counts_cntrls[subgroup][gt_category][(i, b), (b, i)] = 0
    #Example of passing_d:{gt_category:{GT:count}}: {'obs_homo1': {((0, 0), (0, 0)): 0}, 'obs_het': {((0, 1), (1, 0)): 0, ((0, 2), (2, 0)): 0, ((0, 3), (3, 0)): 0, ((1, 2), (2, 1)): 0, (1)): 0, ((2, 3), (3, 2)): 0}, 
#'obs_homo2': {((1, 1), (1, 1)): 0, ((2, 2), (2, 2)): 0, ((3, 3), (3, 3)): 0}}    
    for k,sm in samples.items():
        subgroup = mi.sa.sa_collection[k].get_subgroup()
        if (None in sm['GT']):
            missing += 1
            sm['GT'] = (None, None)
            tallyMissingSample(k, sm)
            continue
        try:
            if (sm['DP'] < cfg.MINDP or sm['GQ'] < cfg.MINGQ):
                if sm['GT'][0] == sm['GT'][1]: #Homozygous Sample (REF or ALT)  
                    if sm['GT'] == (0,0):
                        failing_d['obs_homo1'][((0, 0), (0, 0))]  += 1
                        mi.sa.sa_collection[k].tallySA['failing_obs_homo1'] += 1    # Add Failing samples to SA
                    else: #Homozygous ALT
                        failing_d['obs_homo2'][(sm['GT'], (sm['GT'][1],sm['GT'][0]))] +=1
                        mi.sa.sa_collection[k].tallySA['failing_obs_homo2'] += 1
                else: #Heterozygous sample
                    failing_d['obs_het'][(sm['GT'], (sm['GT'][1],sm['GT'][0]))] +=1
                    if 0 in sm['GT']: # Treat Alternate alleles in a Het sample differently for Zhet calculations
                        mi.sa.sa_collection[k].tallySA['failing_obs_het'] += 1
                    else:
                        mi.sa.sa_collection[k].tallySA['failing_obs_homo2'] += 1

                sm['GT'] = (None, None)
                gt_failed += 1
                tallyFailedSample(k, sm)
                depth_sum += sm['DP']
                continue
        except TypeError:
            for classification in failing_d.keys():
                    for key, values in failing_d[classification].items():
                        if sm['GT'] in key:
                            failing_d[classification][key] += 1
                            mi.sa.sa_collection[k].tallySA['failing_'+ classification] += 1
            sm['GT'] == (None, None)
            tallyFailedSample(k, sm)
            gt_failed += 1
            continue
        except:
            raise
        depth_sum += sm['DP']
        #Increase GT Counts if GT found in sample
        if sm['GT'][0] == sm['GT'][1]: #Homozygous Sample (REF or ALT)
            if sm['GT'] == (0,0):
                passing_d['obs_homo1'][((0, 0), (0, 0))] +=1
                mi.sa.sa_collection[k].tallySA['passing_obs_homo1'] += 1
            else: #Homozygous ALT
                passing_d['obs_homo2'][(sm['GT'], (sm['GT'][1],sm['GT'][0]))] +=1
                mi.sa.sa_collection[k].tallySA['passing_obs_homo2'] += 1
        else: #Heterozygous sample
            passing_d['obs_het'][(sm['GT'], (sm['GT'][1],sm['GT'][0]))] +=1
            if 0 in sm['GT']: # Treat Alternate alleles in a Het sample differently for Zhet calculations
                mi.sa.sa_collection[k].tallySA['passing_obs_het'] += 1
            else:
                mi.sa.sa_collection[k].tallySA['passing_obs_homo2'] += 1
            for allele in sm['GT']:
                abhet_AD_list[allele] += sm['AD'][allele]
            abhet_DP_list[sm['GT'][0]] += ( sm['AD'][sm['GT'][0]] + sm['AD'][sm['GT'][1]] )
            abhet_DP_list[sm['GT'][1]] += ( sm['AD'][sm['GT'][1]] + sm['AD'][sm['GT'][0]] )

        tallyPassingSample(k,sm)


    clean_passing_d = copy.deepcopy(passing_d)
    return [passing_d,failing_d,missing,gt_failed,clean_passing_d,depth_sum,abhet_AD_list,abhet_DP_list]

def count_gt_chrx(male_samples,female_samples,rec_details):
    N = len(rec_details['alt'])
    allele_list = [n for n in range(0,N+1)]
    #allele_list = [0, 1, 2, 3]
    missing = 0
    depth_sum = 0
    gt_failed = 0
    abhet_AD_list = [0 for i in allele_list]
    abhet_DP_list = copy.deepcopy(abhet_AD_list)
    mi.sa.clear_mpairs()
    passing_d_male = {'obs_homo1':{},'obs_het':{},'obs_homo2':{}}
    failing_d_male = {'obs_homo1':{},'obs_het':{},'obs_homo2':{}}
    passing_d_female = {'obs_homo1':{},'obs_het':{},'obs_homo2':{}}
    failing_d_female = {'obs_homo1':{},'obs_het':{},'obs_homo2':{}}
    clean_d = {'male':0,'female':0}
    
    for i in allele_list:    #Create Dictionary with keys being all possible Genotypes
            for b in allele_list[i:]:
                if (i,b) == (0,0) and (b,i) == (0,0):
                    passing_d_male['obs_homo1'][(i,b),(b,i)] = 0
                    failing_d_male['obs_homo1'][(i,b),(b,i)] = 0
                    passing_d_female['obs_homo1'][(i,b),(b,i)] = 0
                    failing_d_female['obs_homo1'][(i,b),(b,i)] = 0
                elif (i,b) != (b,i):
                    passing_d_male['obs_het'][(i,b),(b,i)] = 0
                    failing_d_male['obs_het'][(i,b),(b,i)] = 0
                    passing_d_female['obs_het'][(i,b),(b,i)] = 0
                    failing_d_female['obs_het'][(i,b),(b,i)] = 0
                elif (i,b) == (b,i):
                    passing_d_male['obs_homo2'][(i,b),(b,i)] = 0
                    failing_d_male['obs_homo2'][(i,b),(b,i)] = 0
                    passing_d_female['obs_homo2'][(i,b),(b,i)] = 0
                    failing_d_female['obs_homo2'][(i,b),(b,i)] = 0

    #Example of passing_d:{classification:{GT:count}}: {'obs_homo1': {((0, 0), (0, 0)): 0}, 'obs_het': {((0, 1), (1, 0)): 0, ((0, 2), (2, 0)): 0, ((0, 3), (3, 0)): 0, ((1, 2), (2, 1)): 0, (1)): 0, ((2, 3), (3, 2)): 0}, 
#'obs_homo2': {((1, 1), (1, 1)): 0, ((2, 2), (2, 2)): 0, ((3, 3), (3, 3)): 0}}    
    ref = rec_details['ref']
    alt = rec_details['alt'][0]
    for k,sm in male_samples.items():
        
        subgroup = mi.sa.sa_collection[k].get_subgroup()
        if (None in sm['GT']):
            missing += 1
            sm['GT'] = (None, None)
            tallyMissingSample(k, sm)
            continue
        try:
            if (sm['DP'] < cfg.MINDP or sm['GQ'] < cfg.MINGQ):
                for classification in failing_d_male.keys():
                    for key, values in failing_d_male[classification].items():
                        if sm['GT'] in key:
                            failing_d_male[classification][key] += 1
                            if classification == "obs_homo1":
                                mi.sa.sa_collection[k].tallySA['failing_'+ classification] += 1    # Add Failing samples to SA
                            else:
                                if 0 in sm['GT']:
                                    mi.sa.sa_collection[k].tallySA['failing_obs_het'] += 1
                                else:
                                    mi.sa.sa_collection[k].tallySA['failing_obs_homo2'] += 1
                sm['GT'] = (None, None)
                tallyFailedSample(k, sm)
                gt_failed += 1
                depth_sum += sm['DP']
                continue
        except TypeError:
            for classification in failing_d_male.keys():
                    for key, values in failing_d_male[classification].items():
                        if sm['GT'] in key:
                            failing_d_male[classification][key] += 1
                            mi.sa.sa_collection[k].tallySA['failing_'+ classification] += 1
            sm['GT'] == (None,None)
            tallyFailedSample(k, sm)
            depth_sum += sm['DP']
            gt_failed += 1
            continue
        except:
            raise
        depth_sum += sm['DP']
        #Increase GT Counts if GT found in sample
        for classification in passing_d_male.keys():
            for key, values in passing_d_male[classification].items():
                if sm['GT'] in key:
                    passing_d_male[classification][key] += 1
                    if classification == 'obs_homo1':
                        mi.sa.sa_collection[k].tallySA['passing_'+ classification] += 1
                        if mi.sa.sa_collection[k].is_control():
                            pass

                    elif classification == 'obs_het':
                        gt_failed +=1
                        mi.sa.sa_collection[k].tallySA['passing_'+ classification] += 1
                        if 0 in sm['GT']:
                            mi.sa.sa_collection[k].tallySA['failing_obs_het'] += 1
                        else:
                            mi.sa.sa_collection[k].tallySA['failing_obs_homo2'] += 1

                    elif classification == 'obs_homo2':
                        tallyTiTv(k, sm, ref, alt, None)
                        mi.sa.sa_collection[k].tallySA['passing_obs_homo2'] += 1
                    else:
                        raise TypeError("Weird GT")

        tallyPassingSample(k,sm)
    for k,sm in female_samples.items():
        
        subgroup = mi.sa.sa_collection[k].get_subgroup()
        if (None in sm['GT']):
            missing += 1
            sm['GT'] = (None, None)
            tallyMissingSample(k, sm)
            continue
        try:
            if (sm['DP'] < cfg.MINDP or sm['GQ'] < cfg.MINGQ):
                for classification in failing_d_female.keys():
                    for key, values in failing_d_female[classification].items():
                        if sm['GT'] in key:
                            failing_d_female[classification][key] += 1
                            if classification == "obs_homo1":
                                mi.sa.sa_collection[k].tallySA['failing_'+ classification] += 1    # Add Failing samples to SA
                            else:
                                if 0 in sm['GT']:
                                    mi.sa.sa_collection[k].tallySA['failing_obs_het'] += 1
                                else:
                                    mi.sa.sa_collection[k].tallySA['failing_obs_homo2'] += 1
                sm['GT'] = (None, None)
                tallyFailedSample(k, sm)
                gt_failed += 1
                depth_sum += sm['DP']
                continue
        except TypeError:
            for classification in failing_d_female.keys():
                    for key, values in failing_d_female[classification].items():
                        if sm['GT'] in key:
                            failing_d_female[classification][key] += 1
                            mi.sa.sa_collection[k].tallySA['failing_'+ classification] += 1

            sm['GT'] == (None,None)
            tallyFailedSample(k, sm)
            gt_failed += 1
            continue
        except:
            raise
        depth_sum += sm['DP']
        #Increase GT Counts if GT found in sample
        for classification in passing_d_female.keys():
            for key, values in passing_d_female[classification].items():
                if sm['GT'] in key:
                    passing_d_female[classification][key] += 1
                    if classification == 'obs_homo1':
                        mi.sa.sa_collection[k].tallySA['passing_'+ classification] += 1

                    elif classification == 'obs_het':
                        tallyTiTv(k, sm, ref, alt, None)
                        if 0 in sm['GT']:
                            mi.sa.sa_collection[k].tallySA['passing_obs_het'] += 1
                        else:
                            mi.sa.sa_collection[k].tallySA['passing_obs_homo2'] += 1

                    elif classification == 'obs_homo2':
                        tallyTiTv(k, sm, ref, alt, None)
                        mi.sa.sa_collection[k].tallySA['passing_obs_homo2'] += 1
                    else:
                        raise TypeError("Weird GT")
        tallyPassingSample(k,sm)
        for het_gt in passing_d_male['obs_het']:
            if sm['GT'] in het_gt:
                for allele in sm['GT']:   # For ABHET Calculations:
                    abhet_AD_list[allele] += sm['AD'][allele]
                abhet_DP_list[sm['GT'][0]] += ( sm['AD'][sm['GT'][0]] + sm['AD'][sm['GT'][1]] )
                abhet_DP_list[sm['GT'][1]] += ( sm['AD'][sm['GT'][1]] + sm['AD'][sm['GT'][0]] )
        for het_gt in passing_d_female['obs_het']:
            if sm['GT'] in het_gt:
                for allele in sm['GT']:   # For ABHET Calculations:
                    abhet_AD_list[allele] += sm['AD'][allele]
                abhet_DP_list[sm['GT'][0]] += ( sm['AD'][sm['GT'][0]] + sm['AD'][sm['GT'][1]] )
                abhet_DP_list[sm['GT'][1]] += ( sm['AD'][sm['GT'][1]] + sm['AD'][sm['GT'][0]] )

    clean_d['male'],clean_d['female'] = copy.deepcopy(passing_d_male), copy.deepcopy(passing_d_female)

    return [passing_d_male,failing_d_male,passing_d_female,failing_d_female,missing,gt_failed,clean_d,depth_sum,abhet_AD_list,abhet_DP_list]


def increment_subgroup(k, idx, subgroup_counts, subgroup_counts_cntrls):
    subgroup_counts[ mi.sa.sa_collection[k].get_subgroup() ][idx] += 1
    
    if mi.sa.sa_collection[k].is_control():
        subgroup_counts_cntrls[ mi.sa.sa_collection[k].get_subgroup() ][idx] += 1

    return subgroup_counts, subgroup_counts_cntrls

def increment_subgroup_multiallelic(k, classification, GT, subgroup_counts, subgroup_counts_cntrls):
    subgroup_counts[ mi.sa.sa_collection[k].get_subgroup() ][classification][GT] += 1
    
    if mi.sa.sa_collection[k].is_control():
        subgroup_counts_cntrls[ mi.sa.sa_collection[k].get_subgroup() ][classification][GT] += 1
    return subgroup_counts, subgroup_counts_cntrls

def tallyMissingSample(k, sm):
    mi.sa.tally(k, sm, 0)

def tallyFailedSample(k, sm):
    mi.sa.tally(k, sm, 1)

def tallyPassingSample(k,sm):
    mi.sa.tally(k, sm, -1)

def tallyIndelSample(k, ref, alt):
    if ref==alt: return
    if len(alt) > 1 or len(ref) > 1:
        mi.sa.sa_collection[k].tallySA['non_missing_indel'] += 1
        return

def tallyTiTv(k,sm, ref, alt, wes_flag):
    mi.sa.tallyTiTv(k, ref, alt, wes_flag)

def is_good_gt(sm):
    """
    """
    try:
        if sm['GT'] == (None, None):
            return False
        elif (sm['DP'] < cfg.MINDP
            or sm['GQ'] < cfg.MINGQ
            ):
            return False
    except TypeError:  #TypeError: unorderable types: NoneType() < int() (missing DP)
            return False
    except:
        raise

    return True

# deprecated
def count_gt1(samples, samples_list):
    """
    count_gt - sample genotype {0/0, 0/1, 1/1}
               Apply Genotype-level QC: DP<10, GQ<20; set to ./.
    """
    global MINDP, MINGQ
    obs_hts = 0  # zHet01
    obs_hom1 = 0 # zHet00
    obs_hom2 = 0 # zHet11
    missing = 0
    depth_sum = 0
    ct = 0

    for k,sm in samples.items():
        if k not in samples_list: continue

        try:
            if (sm['DP'] < MINDP
                or sm['GQ'] < MINGQ
                ):
                sm['GT'] = (None, None)
                ct += 1
        except TypeError:  #TypeError: unorderable types: NoneType() < int() (missing DP)
            sm['GT'] = (None, None)
            ct += 1
        except:
            raise

        if (sm['GT'] == (None, None)
            or
            sm['DP'] == None):
            missing += 1
            continue

        depth_sum += sm['DP']

        if (sm['GT'] == (0,1)
        or sm['GT'] == (1,0)):
            obs_hts += 1
        elif sm['GT'] == (0,0):
            obs_hom1 += 1
        elif sm['GT'] == (1,1):
            obs_hom2 += 1
        #else: not counted

    return [obs_hom1, obs_hts, obs_hom2, missing, depth_sum]


