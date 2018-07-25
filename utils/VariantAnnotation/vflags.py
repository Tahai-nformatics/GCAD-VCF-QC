#!/bin/env python3

import config as cfg

from utils.stats.statistical import calc_ExcessHet, calc_pHWE
from utils.stats.count_gt import count_gt, count_gt1

def calcVA(snp_samples, snp_record_filter):
    """
    calcVA - get Variant Annotation; VFLAGS and ABHet
    Variant-level QC
    VFLAG 1: Does variant PASS according to GATK, No=1 (fail)
    VFLAG 2: After genotype-level QC, no?
    VFLAG 3: Monomorphic, yes?
    VFLAG 4: Call Rate <80%, yes?
    VFLAG 5: Mean Depth >500, yes?
    VFLAG 6: Departure from Expected Genotype Distribution: Family data -> Excess Heterozygosity; Unrelated -> Hardy-Weinberg equilibrium if MAF>0.01
    VFLAG 7? is multiallelic==1 
    VFLAG 11? WES
    VFLAG 12? ABHet outside limits
    VFLAG 0: With none of the above
    """
    #global isFam, minTranche, miss_rate, max_dp, hetz_lim1, hetz_lim2, hwe_pval, hwe_maf

    vf = []
    pass_cnt = [0,0,0]
    fail_cnt = [0,0,0]

    # VFLAG 1
    if 'PASS' in snp_record_filter:
        pass_snv = 1
        mp_score = 1
        badcall  = 0
    else:
        for k in snp_record_filter.keys():
            if k.startswith('VQSRTrancheSNP'): # VQSRTrancheSNP99.80to99.90
                low, high = k.replace('VQSRTrancheSNP','').split('to')
                if float(low) >= cfg.minTranche:
                    vf.append(1)
                    pass_snv = 0
                    mp_score = 0
                    badcall  = 1
                else:
                    pass_snv = 1
                    mp_score = 1
                    badcall  = 0

    [obs_hom1, obs_hets, obs_hom2, missing, gt_failed, depth_sum, failed, het_ad, het_dp] = count_gt(snp_samples)
    total = obs_hom1 + obs_hets + obs_hom2 + missing + gt_failed
    non_missing = obs_hom1 + obs_hets + obs_hom2

    # VFLAG 2
    if (missing + gt_failed) == total:
        vf.append(2)

    # VFLAG 3
    if obs_hets == 0:
        if obs_hom2 ==0 or obs_hom1 == 0:
            vf.append(3)

    # VFLAG 4
    callrate = 1 - (missing + gt_failed) / total
    if callrate <= (1 - cfg.miss_rate):
        vf.append(4)
        badcall = 1

    # VFLAG 5
    if non_missing > 0:
        if (depth_sum / non_missing) > cfg.max_dp:
            vf.append(5)

    # VFLAG 6
    maf = 0
    if  non_missing > 0:
        maf = (obs_hets + (2 * obs_hom2)) / (2*(non_missing))
        #if (maf > 0.5):
        #        maf = 1 - maf

    if cfg.isFam:
        z_het, hetz_maf = calc_ExcessHet(obs_hom1, obs_hets, obs_hom2)
        if z_het == '.': z_het = 0

        if   ((maf <  0.2  or maf  > 0.8) and (abs(float(z_het)) > cfg.hetz_lim1)):
            vf.append(6)
        elif ((maf >= 0.2 and maf <= 0.8) and (abs(float(z_het)) > cfg.hetz_lim2)):
            vf.append(6)

    else:
        z_het = '.'

        if non_missing > 0:
            if (maf > 0.5):
                maf = 1 - maf

            # Calc Hardy-Weinberg equilibrium if MAF>0.01
            if maf > hwe_maf:
                z_het = calc_pHWE(obs_hom1, obs_hets, obs_hom2)

        if((z_het >= 1) or (z_het < hwe_pval)):
            vf.append(6)

    # VFLAG 0
    # Presense of VFLAGs counts as failing GTs
    if len(vf) < 1:
        vf.append(0)

        # Set passing
        pass_cnt = [obs_hom1, obs_hets, obs_hom2]
        fail_cnt = failed
    else:
        fail_cnt = [obs_hom1 + failed[0], obs_hets + failed[1], obs_hom2 + failed[2]]

    clean_obs = [obs_hom1, obs_hets, obs_hom2]

    # AB Het
    if het_dp > 0:
        ab_het = "{0:.4f}".format(het_ad / het_dp)
    else:
        ab_het = "0.0000"

    return [vf, ab_het, pass_cnt, fail_cnt, missing, gt_failed, depth_sum, clean_obs]


def calcVFlags1(snp_samples, snp_record_filter, samples_list):
    """
    Variant-level QC
    VFLAG 1: Does variant PASS according to GATK, No=1 (fail)
    VFLAG 2: After genotype-level QC, no?
    VFLAG 3: Monomorphic, yes?
    VFLAG 4: Call Rate <80%, yes?
    VFLAG 5: Mean Depth >500, yes?
    VFLAG 6: Departure from Expected Genotype Distribution: Family data -> Excess Heterozygosity; Unrelated -> Hardy-Weinberg equilibrium if MAF>0.01
    VFLAG 7? is multiallelic==1 
    VFLAG 11? WES
    VFLAG 12? ABHet limit
    """
    global isFam, minTranche, miss_rate, max_dp, hetz_lim1, hetz_lim2, hwe_pval, hwe_maf

    vf = []

    # VFLAG 1
    if 'PASS' in snp_record_filter:
        pass_snv = 1
        mp_score = 1
        badcall  = 0
    else:
        for k in snp_record_filter.keys():
            if k.startswith('VQSRTrancheSNP'): # VQSRTrancheSNP99.80to99.90
                low, high = k.replace('VQSRTrancheSNP','').split('to')
                if float(low) >= minTranche:
                    vf.append(1)
                    pass_snv = 0
                    mp_score = 0
                    badcall  = 1
                else:
                    pass_snv = 1
                    mp_score = 1
                    badcall  = 0
    # VFLAG 2
    [obs_hom1, obs_hets, obs_hom2, missing, depth_sum] = count_gt1(snp_samples, samples_list)
    total = obs_hom1 + obs_hets + obs_hom2 + missing
    if missing == total:
        vf.append(2)

    # VFLAG 3
    if obs_hets == 0:
        if obs_hom2 ==0 or obs_hom1 == 0:
            vf.append(3)

    # VFLAG 4
    callrate = 1 - missing / total
    if callrate <= miss_rate:
        vf.append(4)
        badcall = 1

    # VFLAG 5
    if depth_sum > max_dp:
        vf.append(5)

    # VFLAG 6
    if obs_hom1 + obs_hets + obs_hom2 > 0:
        maf = (obs_hets + (2 * obs_hom2)) / (2*(obs_hom1 + obs_hets + obs_hom2))
        #if (maf > 0.5):
        #        maf = 1 - maf
    else:
        maf = 0

    if isFam:
        z_het, hetz_maf = calc_ExcessHet(obs_hom1, obs_hets, obs_hom2)
        if z_het == '.': z_het = 0

        if   ((maf <  0.2  or maf  > 0.8) and (abs(float(z_het)) > hetz_lim1)):
            vf.append(6)
        elif ((maf >= 0.2 and maf <= 0.8) and (abs(float(z_het)) > hetz_lim2)):
            vf.append(6)

    else:
        z_het = '.'

        if (2 * (obs_hets + obs_hom1 + obs_hom2)) > 0:
            if (maf > 0.5):
                maf = 1 - maf

            # Calc Hardy-Weinberg equilibrium if MAF>0.01
            if maf > hwe_maf:
                z_het = calc_pHWE(obs_hom1, obs_hets, obs_hom2)

        if((z_het >= 1) or (z_het < hwe_pval)):
            vf.append(6)

    return vf


def applyGenotypeQC(samples):
    """
    applyGenotypeQC -
    """
    global MINDP, MINGQ
    ct = 0
    for k,sm in samples.items():

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

    return ct


