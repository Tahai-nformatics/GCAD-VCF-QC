#!/usr/bin/env python3
#-*- coding: utf-8 -*-

from argparse import ArgumentParser
from pysam import VariantFile
from utils.stats.statistical import calc_ExcessHet, calc_pHWE
#import configparser
import csv
from collections import namedtuple

MINDP=10
MINGQ=20
hwe_maf= 0.01
minTranche = 99.7
miss_rate = 0.2
max_dp = 500
isFam = 1
hetz_lim1 = 5
hetz_lim2 = 6.101825
hwe_pval = 5e-06

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

def count_gt(samples):
    """
    count_gt - sample genotype {0/0, 0/1, 1/1}
    """
    obs_hts = 0  # zHet01
    obs_hom1 = 0 # zHet00
    obs_hom2 = 0 # zHet11
    missing = 0
    depth_sum = 0

    for k,sm in samples.items():

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

def extractSubSets(fam):
    """
    @return Samples Dict per Subset
    """
    SampleFamDetail = namedtuple('SampleFamDetail',['SampID','FID','SubjID','FA','MO','SEX','AFF','AD','AGE','ADSPWGS','Subset','Subgroup','Race_Ethnicity','SeqCtr','ExcludeFromZHet'])

    samples = dict()
    for sm in map(SampleFamDetail._make, csv.reader(open(fam, 'r'),delimiter='\t')):
        if sm.Subset in samples:
            samples[sm.Subset].append(sm.SampID)
        else:
            samples[sm.Subset]=[sm.SampID]
    return samples

def calcVFlags(snp_record, gt_failed):
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
    if 'PASS' in snp_record.filter:
        pass_snv = 1
        mp_score = 1
        badcall  = 0
    else:
        for k in snp_record.filter.keys():
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
    [obs_hom1, obs_hets, obs_hom2, missing, depth_sum] = count_gt(snp_record.samples)
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

def main():
    argparser = ArgumentParser()
    grp_file_paths = argparser.add_argument_group(title='File paths')
    grp_file_paths.add_argument('--vcf', type=str, help='input VCF file', required=True)
    grp_file_paths.add_argument('--fam', type=str, help='fam file (with headers)', required=True)
    #grp_file_paths.add_argument('--outfile', type=str, help='filepath for output file', required=True)

    grp_settings = argparser.add_argument_group(title='Settings')
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

    MINDP = args.min_dp
    MINGQ = args.min_gq
    minTranche = args.minTranche
    miss_rate = args.miss_rate
    max_dp = args.max_dp
    isFam = args.isFam
    hetz_lim1 = args.hetz_lim1
    hetz_lim2 = args.hetz_lim2
    hwe_pval = args.hwe_pval
    hwe_maf = args.hwe_maf

    samples = extractSubSets(args.fam)
    #{k:rec.samples[k] for k in samples['sub1']}
    #{k:rec.samples[k] for k in samples['sub1'] if k in rec.samples}
    #{key: d[key] for key in d.viewkeys() & l}

    vcf_in = VariantFile(args.vcf)
    #vcf_out = VariantFile(args.outfile, 'w', header=vcf_in.header, threads=4)



    for rec in vcf_in.fetch():
        #vcf_out.write(rec)

        # Apple Genotype-level QC: DP<10, GQ<20; set to ./.
        gt_failed = applyGenotypeQC(rec.samples)

        vf = calcVFlags(rec, gt_failed)



if __name__ == "__main__":
    main()
