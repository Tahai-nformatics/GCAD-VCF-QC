#!/usr/bin/env python3

import config as cfg
import csv
#from utils.stats.statistical import calc_ExcessHet, calc_pHWE
from utils.stats.count_gt import count_gt

# global dictionary containing target intervals used in WES QC
targets = dict()

def calcVA(snp_samples, rec_details, subset):
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
    VFLAG 11: WES; Does  this  variant  fall  within  the  provided  target  capture  regions, no?
    VFLAG 12? ABHet outside limits
    VFLAG 0: With none of the above
    """

    vf = []
    pass_cnt = [0,0,0]
    fail_cnt = [0,0,0]
    snp_record_filter = rec_details['filter']

    # VFLAG 1
    if 'PASS' in snp_record_filter:
        pass_snv = 1
    else:
        for k in snp_record_filter.keys():
            if k.startswith('VQSRTranche'): # VQSRTrancheSNP99.80to99.90; VQSRTrancheINDEL
                low = None
                if k.startswith('VQSRTrancheSNP'):
                    low, high = k.replace('VQSRTrancheSNP','').split('to')
                elif k.startswith('VQSRTrancheINDEL'):
                    low, high = k.replace('VQSRTrancheINDEL','').split('to')

                if low is not None:
                    if float(low) >= cfg.minTranche:
                        vf.append(1)
                        pass_snv = 0
                    else:
                        pass_snv = 1

    [obs_hom1, obs_hets, obs_hom2, missing, gt_failed, depth_sum, failed, het_ad, het_dp, subg, subg_c] = count_gt(snp_samples, rec_details)
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

    # VFLAG 5
    if non_missing > 0:
        if (depth_sum / non_missing) > cfg.max_dp:
            vf.append(5)

    # VFLAG 6
    #maf = 0
    #if  non_missing > 0:
        #maf = (obs_hets + (2 * obs_hom2)) / (2*(non_missing))
        ##if (maf > 0.5):
        ##        maf = 1 - maf

    #if cfg.isFam:
        #z_het, hetz_maf = calc_ExcessHet(obs_hom1, obs_hets, obs_hom2)
        #if z_het == '.': z_het = 0

        #if   ((maf <  0.2  or maf  > 0.8) and (abs(float(z_het)) > cfg.hetz_lim1)):
            #vf.append(6)
        #elif ((maf >= 0.2 and maf <= 0.8) and (abs(float(z_het)) > cfg.hetz_lim2)):
            #vf.append(6)

    #else:
        #z_het = '.'

        #if non_missing > 0:
            #if (maf > 0.5):
                #maf = 1 - maf

            ## Calc Hardy-Weinberg equilibrium if MAF>0.01
            #if maf > hwe_maf:
                #z_het = calc_pHWE(obs_hom1, obs_hets, obs_hom2)

        #if((z_het >= 1) or (z_het < hwe_pval)):
            #vf.append(6)

    # VFLAG 7
    #if len(rec_details['alt']) > 1:
    #    vf.append(7)

    # VFLAG 11
    if targets and subset in targets[rec_details['chr']]:
        variant_bin = reg2bin(rec_details['pos'])

        # simple case: bin is not in targets dict
        if variant_bin not in targets[rec_details['chr']][subset]:
            vf.append(11)
        else:
            in_region = False
            is_insertion = len(rec_details['alt'][0]) - 1
            is_deletion = len(rec_details['ref']) - 1

            for interval in targets[rec_details['chr']][subset][variant_bin]:

                if is_deletion:

                    lower_bound = rec_details['pos']
                    upper_bound = rec_details['pos'] + is_deletion

                    # signal if the indel is contained within the target
                    if interval[0] <= lower_bound <= interval[1]:
                        in_region = True
                        break

                    if interval[0] <= upper_bound <= interval[1]:
                        in_region = True
                        break

                    # also signal if the target is within the interval
                    if lower_bound <= interval[0] <= upper_bound:
                        in_region = True
                        break

                    if lower_bound <= interval[1] <= upper_bound:
                        in_region = True
                        break

                else:
                    if interval[0] <= rec_details['pos'] <= interval[1]:
                        in_region = True
                        break

            if not in_region:
                vf.append(11)

    # VFLAG 0
    # Presence of VFLAGs counts as failing GTs
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

    return [vf, ab_het, pass_cnt, fail_cnt, missing, gt_failed, depth_sum, clean_obs, subg, subg_c]


def read_target_files(target_list, chr):
    """
    read_target_files - read-in BED files for WES VFLAG 11
    :param targets: list of BED files
    :param chr: optional chromosome region restriction
    """

    bin_set = set()

    for trgt_str in target_list:
        trgt = trgt_str.split(':')

        if len(trgt) > 1:
            subset = trgt[1]

            if chr not in targets:
                targets[chr] = {subset: dict()}

        else:
            raise ValueError("WES target file missing subset assignment")

        with open(trgt[0]) as bed_file:

            for bed_line in bed_file:
                row = bed_line.split()
                f_chr = row[0]
                f_start = int(row[1]) - cfg.flank_size
                f_end = int(row[2]) + cfg.flank_size
                if chr:
                    if f_chr != chr: continue
                else:
                    if f_chr not in targets:
                        targets[f_chr] = {subset: dict()}

                bin_start = reg2bin(f_start)
                bin_end = reg2bin(f_end)

                bin_set.clear()

                if bin_start == bin_end:
                     bin_set.add(bin_start)
                else:
                   bin_set.add(bin_start)
                   bin_set.add(bin_end)
                   for i in range(f_start + 1, f_end):
                      bin_set.add ( reg2bin(i) )


                if subset not in targets[f_chr]:
                    targets[f_chr][subset] = dict()

                for region_bin in bin_set:
                   if region_bin in targets[f_chr][subset]:
                       targets[f_chr][subset][region_bin].append([f_start, f_end])
                   else:
                       targets[f_chr][subset][region_bin] = list([[f_start, f_end]])

    return


def reg2bin(beg):
    """
    reg2bin - convert region to bin, adapted from genomic interval conversion to bin position
              Based off the algorithm presented in:
              https://samtools.github.io/hts-specs/SAMv1.pdf

    :param beg: region start
    :return: bin
    """
    #end = beg
    #if beg >> 14 == end >> 14: return int(((1 << 15)-1) / 7 + (beg >> 14))
    #if beg >> 17 == end >> 17: return int(((1 << 12)-1) / 7 + (beg >> 17))
    #if beg >> 20 == end >> 20: return int(((1 << 9)-1) / 7 + (beg >> 20))
    #if beg >> 23 == end >> 23: return int(((1 << 6)-1) / 7 + (beg >> 23))
    #if beg >> 26 == end >> 26: return int(((1 << 3)-1) / 7 + (beg >> 26))
    return int(((1 << 15)-1) / 7 + (beg >> 14))

