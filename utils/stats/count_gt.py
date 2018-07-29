#!/usr/bin/env python3

import config as cfg
import utils.SampleAnnotation.sample_annotation as mi

def count_gt(samples,rec):
    """
    count_gt - sample genotype {0/0, 0/1, 1/1}
               Apply Genotype-level QC: DP<10, GQ<20; set to ./.
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
    singleton = ''
    doubletons = []
    ref = rec.ref
    alt = rec.alts[0]

    mi.sa.clear_mpairs()

    for k,sm in samples.items():


        if (sm['GT'] == (None, None)
            ):
            missing += 1
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

                #sm['GT'] += (0)
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

            #sm['GT'] += (0)
            tallyFailed(k, sm)
            gt_failed += 1
            #raise TypeError("Weird {},{}, {}, {}".format(k,str(sm['GT']), str(sm['DP']), str(sm['AD']) ))
            continue
        except:
            raise



        depth_sum += sm['DP']

        if (sm['GT'] in {(0,1), (1,0)} ):
            obs_hts += 1
            het_ad += sm['AD'][0]
            het_dp += sm['DP']
        elif sm['GT'] == (0,0):
            obs_hom1 += 1
        elif sm['GT'] == (1,1):
            obs_hom2 += 1
        else:
            raise TypeError("Weird GT")

        tallyPassing(k, sm, ref, alt)

        # Record possible singleton, doubletons
        if (obs_hts + obs_hom2) == 1:
            singleton = k
        elif (obs_hts + obs_hom2) == 2:
            doubletons = [singleton, k]
            #singleton = None


    if (obs_hts + obs_hom2) == 1:
        mi.sa.sa_collection[singleton].tallySA['singleton'] += 1
    if (obs_hts + obs_hom2) == 2:
        for indiv in doubletons:
            mi.sa.sa_collection[indiv].tallySA['doubleton'] += 1

    return [obs_hom1, obs_hts, obs_hom2, missing, gt_failed, depth_sum, failed, het_ad, het_dp]

def tallyMissing(k, sm):
    mi.sa.tally(k, sm, 0)

def tallyFailed(k, sm):
    mi.sa.tally(k, sm, 1)

def tallyPassing(k,sm, ref, alt):
    mi.sa.tally(k, sm, -1)
    mi.sa.tallyTiTv(k, ref, alt)

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


