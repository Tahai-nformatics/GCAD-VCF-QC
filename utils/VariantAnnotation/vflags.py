#!/usr/bin/env python3

import config as cfg
import csv
import copy
from scipy.stats import binom
import scipy.stats as stats
from utils.stats.count_gt import count_gt, count_gt_multiallelic, count_gt_chrx
# global dictionary containing target intervals used in WES QC
targets = dict()

# global dictionart containing exons
exons = dict()



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

    # VFLAG 11
    in_region = None
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

    in_exon = check_inside_exon(rec_details['pos'], rec_details['chr'])

    # sample-level qc below. Results are required for remaining vflags
    [obs_hom1, obs_hets, obs_hom2, missing, gt_failed, depth_sum, failed, het_ad, het_dp ] = count_gt(snp_samples, rec_details, in_exon)
    total = obs_hom1 + obs_hets + obs_hom2 + missing + gt_failed
    non_missing = obs_hom1 + obs_hets + obs_hom2
    total_genotypes = non_missing + gt_failed
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
    if total_genotypes > 0:
        if (depth_sum / total_genotypes) > cfg.max_dp:
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

    # VFLAG 0
    # Presence of VFLAGs counts as failing GTs
    if len(vf) < 1:
        vf.append(0)

        # Set passing
    pass_cnt = [obs_hom1, obs_hets, obs_hom2]
    fail_cnt = failed
    #else:
    #    fail_cnt = [obs_hom1 + failed[0], obs_hets + failed[1], obs_hom2 + failed[2]]

    #clean_obs = [obs_hom1, obs_hets, obs_hom2]

    # AB Het
    if het_dp > 0:
        ab_het = "{0:.4f}".format(het_ad / het_dp)
    else:
        ab_het = '.'

    return [vf, ab_het, pass_cnt, fail_cnt, missing, gt_failed, depth_sum]


def calcVA_multiallelic(snp_samples,rec_details,subset,vtype):
    """
    """

    vf = []
    ab_het = 0
    total = 0
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

    #Skipping VLAG 11 (WES)

    [passing_d,failing_d,missing,gt_failed,clean_passing_d,depth_sum,abhet_AD_list,abhet_DP_list, allele_count_dict]= count_gt_multiallelic(snp_samples,rec_details,vtype)
    obs_hom1 = sum(list(passing_d['obs_homo1'].values()))
    obs_het = sum(list(passing_d['obs_het'].values()))
    obs_hom2 = sum(list(passing_d['obs_homo2'].values()))
    ab_het = [i for i in range(len(abhet_AD_list))]
    total = obs_hom1 + obs_het + obs_hom2 + missing + gt_failed
    non_missing = obs_hom1+obs_het+obs_hom2
    total_genotypes = non_missing + gt_failed
    N = len(rec_details['alt'])
    allele_list = [n for n in range(0,N+1)]
    maf = []
    sum_clean = sum(list(clean_passing_d['obs_homo1'].values())) + sum(list(clean_passing_d['obs_het'].values())) + sum(list(clean_passing_d['obs_homo2'].values()))
    
    #Calculate MAF before vflag assigned to Passing_d and samples fail
    het_maf_dict = {}
    homo_maf_dict = {}
    ac_ref_het = 0
    temp = sum_clean *2
    if temp >0:
        for allele in allele_list:
            het_maf_dict[allele] = 0
            homo_maf_dict[allele] = 0
            for key in passing_d['obs_het']:
                if allele in key[0]:
                    ac_ref_het += 1
                    het_maf_dict[allele] += passing_d['obs_het'][key]
            for key in passing_d['obs_homo2']:
                if allele in key[0]:
                    homo_maf_dict[allele] += passing_d['obs_homo2'][key]
                else:
                    continue
            if allele == 0:
                maf.append(float(("{0:.6f}".format((het_maf_dict[allele] + (2 * sum(list(passing_d['obs_homo1'].values())))) / temp))))
            elif allele!=0:
                maf.append(("{0:.6f}".format((het_maf_dict[allele] + (2 * homo_maf_dict[allele])) / temp)))
    else: #All samples Missing or failed
        for allele in allele_list:
            maf.append(format(0.0, '.6f'))

    # VFLAG 2
    if (missing + gt_failed) == total:
        vf.append(2)
    # VFLAG 3
    if obs_het ==0:
        if obs_hom1 ==0 or obs_hom2 ==0:
            ct =0
            for i in list(clean_passing_d['obs_homo2'].values()):
                if int(i) > 0:
                    ct +=1
            if ct >=2:
                pass
            else:
                vf.append(3)

    #VFLAG 4 
    callrate = 1 - (missing + gt_failed) / total
    if callrate <= (1-cfg.miss_rate):
        vf.append(4)
    # VFlag 5:
    if total_genotypes > 0:
        if (depth_sum / total_genotypes) > cfg.max_dp:
            vf.append(5)

    #ABHet
    if sum(abhet_DP_list) > 0:
        for allele in ab_het:
            if abhet_DP_list[allele] == 0:
                ab_het[allele] = '.'
            else:
                ab_het[allele] = "{0:.5f}".format(abhet_AD_list[allele] / abhet_DP_list[allele])
                if ab_het[allele] == '0.0000':
                    ab_het[allele] = '.'
    else:
        for allele in ab_het:
            ab_het[allele] = '.'
    

    if len(vf) < 1:
        vf.append(0)

    for item in range(len(ab_het)):
        try:
            ab_het[item] = float(ab_het[item])
        except:
            pass

    return [vf,maf,passing_d,failing_d,missing,gt_failed,clean_passing_d,sum_clean,depth_sum,ab_het, allele_count_dict]

def calcVA_chrx(male_snp_samples,female_snp_samples,rec_details,male_subset,female_subset, chrx_is_multiallelic, vtype):
    vf = []
    ab_het = 0
    total = 0
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
    #Skipping VLAG 11 (WES)
    [passing_d_male,failing_d_male,passing_d_female,failing_d_female,missing,gt_failed,clean_d,depth_sum,abhet_AD_list,abhet_DP_list,allele_count_dict]= count_gt_chrx(male_snp_samples,female_snp_samples,rec_details, chrx_is_multiallelic, vtype)
    ab_het = [i for i in range(len(abhet_AD_list))]
    obs_hom1_male = sum(list(passing_d_male['obs_homo1'].values()))
    obs_hom1_female = sum(list(passing_d_female['obs_homo1'].values()))
    obs_het_male = sum(list(passing_d_male['obs_het'].values()))
    obs_het_female = sum(list(passing_d_female['obs_het'].values()))
    obs_hom2_male = sum(list(passing_d_male['obs_homo2'].values()))
    obs_hom2_female = sum(list(passing_d_female['obs_homo2'].values()))
    allele_list = [n for n in range(0,len(rec_details['alt'])+1)]
    maf_male = copy.deepcopy(clean_d['male']) #maf_male Het GT's will be set to 0, and used in maf calculation
    maf_female = copy.deepcopy(clean_d['female'])
    maf = []
    maf_reference_alleles = sum(list(clean_d['male']['obs_homo1'].values())) + 2*sum(list(clean_d['female']['obs_homo1'].values()))        
    temp = 2 * (sum(list(clean_d['female']['obs_homo1'].values())) + sum(list(clean_d['female']['obs_homo2'].values())) + sum(list(clean_d['female']['obs_het'].values()))) + sum(list(clean_d['male']['obs_homo1'].values())) + sum(list(clean_d['male']['obs_homo2'].values()))
    sum_clean = 0

#Add All Non-Male_Het GT's to sum_Clean
    for key1,key2 in zip(clean_d['male'].keys(),clean_d['female'].keys()): #obs_hom1, obs_het, obs_hom2
        if key1 != 'obs_het':
            for val1,val2 in zip(clean_d['male'][key1].values(), clean_d['female'][key2].values()): #Genotypes 
                sum_clean += val1 + val2                   #add values to sum_clean
        else: #Add female Het GT's
            for value in clean_d['female'][key1].values(): #Genotypes
                sum_clean += value
    
    total = missing + gt_failed + sum_clean
#For MAF calculation, set Male_Passing_Het to 0
    for k,v in maf_male['obs_het'].items():
        maf_male['obs_het'][k] = 0

     #Calculate MAF before vflag assigned to Passing_d and samples fail
    het_maf_dict = {}
    homo_maf_dict = {}
    ac_ref_het = 0

    if temp >0:
        for allele in allele_list:
            het_maf_dict[allele] = 0
            homo_maf_dict[allele] = 0
            for key_male,key_female in zip(maf_male['obs_het'],maf_female['obs_het']): # i.e key_male = ((0, 1), (1, 0)) key_male[0] = (0,1)
                if allele in key_male[0]:
                    het_maf_dict[allele] += maf_male['obs_het'][key_male]
                if allele in key_female[0]:
                    het_maf_dict[allele] += maf_female['obs_het'][key_female]
            for key_male,key_female in zip(maf_male['obs_homo2'],maf_female['obs_homo2']):
                if allele in key_male[0]:
                    homo_maf_dict[allele] += maf_male['obs_homo2'][key_male]
                if allele in key_female[0]:
                    homo_maf_dict[allele] += 2*maf_female['obs_homo2'][key_female] #2*maf_female because 2 alleles for female
                else:
                    continue
            if allele == 0:
                maf.append(float(("{0:.6f}".format((het_maf_dict[allele] + (  maf_reference_alleles)) / temp))))
            elif allele!=0:
                maf.append(("{0:.6f}".format((het_maf_dict[allele] + ( homo_maf_dict[allele])) / temp)))
    else:
        for allele in allele_list:
            maf.append(format(0.0, '.6f'))


    # VFLAG 2
    if (missing + gt_failed) == total:
        vf.append(2)
    
    # VFLAG 3
    total_obs_hom1 = obs_hom1_male + obs_hom1_female
    total_obs_hom2 = obs_hom2_male + obs_hom2_female
    if obs_het_female ==0:
        if total_obs_hom1 ==0 or total_obs_hom2 ==0:
            ct =0
            for i in list(clean_d['male']['obs_homo2'].values()):
                if int(i) > 0:
                    ct +=1
            if ct >=2:
                pass
            else:
                vf.append(3)

    #VFLAG 4
    callrate = 1 - (missing + gt_failed) / total
    if callrate <= (1-cfg.miss_rate):

        vf.append(4)
    # VFlag 5:
    if total > 0:
        if (depth_sum / total) > cfg.max_dp:
            vf.append(5)

    if sum(abhet_DP_list) > 0:
        for item in ab_het:
            if abhet_DP_list[item] == 0:
                ab_het[item] = '.'
            else:
                ab_het[item] = "{0:.4f}".format(abhet_AD_list[item] / abhet_DP_list[item])
                if ab_het[item] == '0.0000':
                    ab_het[item] = '.'
    else:
        ab_het = '.'
    
    #VLFAG 7 if any Male_Het GT is > 6 then set variant to VFLAG 7
    t=0.0001 #Prob of false positive
    e=0.0001 #error rate
    N= passing_d_male['obs_het'].values()
    N = len(male_snp_samples)
    c = stats.binom.ppf((1-t),N,e)
    # qbinom((1-t),N,e) #gives number c such that P(bin>c)<t. We reject a position if #male hets > c (strictly greater than)
    if sum(passing_d_male['obs_het'].values()) > c:
        vf.append(7)

    #Add the Male_Hets which "passed" to failing_d_male['obs_het'] and add that to gt_failed
    if len(vf) < 1:
        vf.append(0)
        for key, values in passing_d_male['obs_het'].items():
            failing_d_male['obs_het'][key] += values

    for item in range(len(ab_het)):
        try:
            ab_het[item] = float(ab_het[item])
        except:
            pass
    
    return [vf,passing_d_male,passing_d_female,failing_d_male,failing_d_female,missing,gt_failed,clean_d,sum_clean,maf, depth_sum,ab_het,allele_count_dict]

def check_inside_exon(pos, contig):
   """
   check_inside_exon - check to see if the variant position is inside an exon in exon_file
   return true if pos is within start-end, false otherwise
   """
   if exons and contig in exons:
        variant_bin = reg2bin(pos)

        # simple case: bin is not in exons dict
        if variant_bin not in exons[contig]:
            return False
        else:

            for interval in exons[contig][variant_bin]:
               if interval[0] <= pos <= interval[1]:
                   return True

   return False

def read_exon_file(exon_file, chr):
   """
   read_exon_file - read-in exon file from BED format. Store as a dictionary of bins as keys
                    and start/end as values
   """

   bin_set = set() # data structure to temporarily hold bin values

   if chr:
     exons[chr] = {}


   with open(exon_file) as bed_file:

       for bed_line in bed_file:
           row = bed_line.split()
           f_chr = row[0]
           f_start = int(row[1])
           f_end = int(row[2])

           if chr:
               if f_chr != chr: continue
           else:
               if f_chr not in exons:
                  exons[f_chr] = {}

           bin_set.clear()

           bin_start = reg2bin(f_start)
           bin_end = reg2bin(f_end)

           if bin_start == bin_end:
               bin_set.add(bin_start)
           else:
               bin_set.add(bin_start)
               bin_set.add(bin_end)
               for i in range(f_start + 1, f_end):
                  bin_set.add ( reg2bin(i) )

           for region_bin in bin_set:
               if region_bin in exons[f_chr]:
                   exons[f_chr][region_bin].append([f_start, f_end])
               else:
                   exons[f_chr][region_bin] = list([[f_start, f_end]])



def read_target_files(target_list, chr):
    """
    read_target_files - read-in BED files for WES VFLAG 11. We build a targets(dict) to store
                        all lines within each BED file. The keys are chromosomes. Each chromosome
                        key has a dict() of subsets listing each bin for associated target region line.
                        So all target region lines are converted to bin and grouped by their bin, and
                        then grouped by their subset.

    :param target_list: list of BED files formatted having a colon ':' between the path and subset name
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
            print("Reading %s" % trgt[0])

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



