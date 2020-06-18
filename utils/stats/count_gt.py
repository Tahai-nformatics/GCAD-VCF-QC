#!/usr/bin/env python3

import config as cfg
import utils.SampleAnnotation.sample_annotation as mi
from collections import OrderedDict
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

    subgroup_counts = OrderedDict({key:[0,0,0] for key in mi.sa.subgroups})
    subgroup_counts_cntrls = copy.deepcopy(subgroup_counts)

    ref = rec_details['ref']
    alt = rec_details['alt'][0]

    mi.sa.clear_mpairs()

    for k,sm in samples.items():


        if (None in sm['GT']
            ):
            missing += 1
            sm['GT'] = (None, None)
            tallyMissing(k, sm)
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
                tallyFailed(k, sm)
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
            tallyFailed(k, sm)
            gt_failed += 1
            #raise TypeError("Weird {},{}, {}, {}".format(k,str(sm['GT']), str(sm['DP']), str(sm['AD']) ))
            continue
        except:
            raise



        depth_sum += sm['DP']
        tallyIndel(k, ref, alt)

        if (sm['GT'] in {(0,1), (1,0)} ):
            obs_hts += 1
            het_ad += sm['AD'][0]
            het_dp += sm['DP']

            subgroup_counts, subgroup_counts_cntrls = increment_subgroup(k, 1, subgroup_counts, subgroup_counts_cntrls)
            tallyTiTv(k, sm, ref, alt, in_region)
        elif sm['GT'] == (0,0):
            obs_hom1 += 1
            subgroup_counts, subgroup_counts_cntrls = increment_subgroup(k, 0, subgroup_counts, subgroup_counts_cntrls)
        elif sm['GT'] == (1,1):
            obs_hom2 += 1
            subgroup_counts, subgroup_counts_cntrls = increment_subgroup(k, 2, subgroup_counts, subgroup_counts_cntrls)
            #if obs_hom2 == 1 and obs_hts == 0:
                #p_dblton = k
            tallyTiTv(k, sm, ref, alt, in_region)
        else:
            raise TypeError("Weird GT")

        tallyPassing(k, sm)


    return [obs_hom1, obs_hts, obs_hom2, missing, gt_failed, depth_sum, failed, het_ad, het_dp, subgroup_counts, subgroup_counts_cntrls]

def increment_subgroup(k, idx, subgroup_counts, subgroup_counts_cntrls):
    subgroup_counts[ mi.sa.sa_collection[k].get_subgroup() ][idx] += 1

    if mi.sa.sa_collection[k].is_control():
        subgroup_counts_cntrls[ mi.sa.sa_collection[k].get_subgroup() ][idx] += 1

    return subgroup_counts, subgroup_counts_cntrls

def tallyMissing(k, sm):
    mi.sa.tally(k, sm, 0)

def tallyFailed(k, sm):
    mi.sa.tally(k, sm, 1)

def tallyPassing(k,sm):
    mi.sa.tally(k, sm, -1)

def tallyIndel(k, ref, alt):
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


