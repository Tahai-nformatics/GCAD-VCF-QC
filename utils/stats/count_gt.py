#!/bin/env python3

import config as cfg

def count_gt(samples):
    """
    count_gt - sample genotype {0/0, 0/1, 1/1}
               Apply Genotype-level QC: DP<10, GQ<20; set to ./.
    """
    #global MINDP, MINGQ
    obs_hts = 0  # zHet01
    obs_hom1 = 0 # zHet00
    obs_hom2 = 0 # zHet11
    missing = 0
    depth_sum = 0
    ct = 0

    for k,sm in samples.items():

        try:
            if (sm['DP'] < cfg.MINDP
                or sm['GQ'] < cfg.MINGQ
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


