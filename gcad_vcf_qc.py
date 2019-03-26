#!/usr/bin/env python3
#-*- coding: utf-8 -*-

from argparse import ArgumentParser
from pysam import VariantFile

import csv
import os
from subprocess import check_output

from collections import namedtuple, OrderedDict, Counter
import time
import cProfile

from utils.VariantAnnotation.vflags import calcVA, read_target_files
from utils.stats.count_gt import is_good_gt
from utils.stats.statistical import calc_ExcessHet, calc_pHWE

import utils.SampleAnnotation.sample_annotation as mi
import config as cfg

import warnings
warnings.simplefilter('always')


def extract_subsets(fam):
    """
    extract_subsets - Creates an OrderedDict() where Subset groupings are the key,
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
        mean_depth = depth_sum / sum_clean if sum_clean else 0

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
            'QUAL':"{0:.2f}".format(rec.qual),'FILTER':",".join(rec.filter.keys()),
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


def write_indiv_summary(prefix):
    """
    """
    outfile = '{}.tsv'.format(prefix)

    with open(outfile, 'w') as csvfile:
        fieldnames = ['SampleID','SEX',
                      'total_nRR','total_nRA','total_nAA','Missing','Set_Missing',
                      'Singleton','Private_Doubleton','Doubleton','HetHom',
                      'Ti','Tv','TiTvRatio','IndMeanDepth',
                      '1P_MI','2P_MI','MI_pairs','Non_Missing_Indels',]
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

            writer.writerow({'SampleID': indiv, 'SEX': val.details_dict.SEX,
                            'total_nRR': val.tallySA[(0,0)],'total_nRA': good_het_gt,'total_nAA': val.tallySA[(1,1)],
                            'Missing': val.tallySA[(None,None)],'Set_Missing': val.tallySA[-9],
                            'Singleton': val.tallySA['singleton'],
                            'Private_Doubleton': val.tallySA['p_dblton'],
                            'Doubleton': val.tallySA['doubleton'],
                            'HetHom':"{0:.5f}".format(het_hom),
                            'Ti': val.tallySA['ti'], 'Tv': val.tallySA['tv'], 'TiTvRatio':"{0:.5f}".format(ti_tv),'IndMeanDepth':"{0:.5f}".format(mean_depth),
                            '1P_MI': val.tallySA['vp1'],'2P_MI':val.tallySA['vp2'],'MI_pairs':val.tallySA['mend_pair'],
                            'Non_Missing_Indels': val.tallySA['non_missing_indel']
                            })


def delete_previous_outputs(out_dir, prefix, subsets):
    """
    """
    mi_file = out_dir + 'summary.mi.tsv'
    for out_file in [mi_file] + [ '{}{}.{}.tsv'.format(out_dir, prefix, x) for x in subsets]:
        if os.path.exists(out_file):
            os.remove(out_file)
    return


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
    grp_settings.add_argument('--min_gq', type=int, help='maximum genotype quality (GQ)', default=20, required=False)
    grp_settings.add_argument('--is_fam', type=int, help='is this family data, 1=yes, 0=no', default=1, required=False)
    grp_settings.add_argument('--min_tranche', type=int, help='Minimum Tranche Score (Filter from VQSR)', default=99.7, required=False)
    grp_settings.add_argument('--miss_rate', type=int, help='max missingness threshold', default=0.2, required=False)
    grp_settings.add_argument('--max_dp', type=int, help='max DP threshold', default=500, required=False)
    grp_settings.add_argument('--hetz_lim1', type=int, help='', default=99999, required=False)
    grp_settings.add_argument('--hetz_lim2', type=int, help='', default=99999, required=False)
    grp_settings.add_argument('--hwe_pval', type=int, help='', default=5e-06, required=False)
    grp_settings.add_argument('--hwe_maf', type=int, help='MAF threshold', default=0.01, required=False)

    wes_settings = argparser.add_argument_group(title='WES Settings')
    wes_settings.add_argument('--wes_target', type=str, help='WES Target BED file, e.g. filename:subset ', nargs='*', required=False)
    wes_settings.add_argument('--flank_size', type=int, help='size of target interval expansion in bp', default=7, required=False)

    args, xtra = argparser.parse_known_args()
    print(args)
    print(xtra)

    ct = 0
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

    # organize new vcf_out header
    vcf_out_hdr = vcf_in.header

    for k in samplesDict.keys():
        vcf_out_hdr.add_meta('INFO', items=[('ID', 'VFLAGS_' + k), ('Number','.'), ('Type', 'String'), ('Description','Pipeline-specific QC variant flags')])
        vcf_out_hdr.add_meta('INFO', items=[('ID', 'ABHet_' + k), ('Number',1), ('Type', 'Float'), ('Description','Allelic Read Ratio')])

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
        print("No output VCF file.")
        print("{}".format(args))
    else:
        baseStr = os.path.basename(args.vcf.replace('.g.vcf','').replace('.vcf','').rpartition('.')[0])
        vcf_out_filename = "{}{}".format(args.out_dir, 'flagged.' + baseStr + regionStr + '.g.vcf.gz')
        vcf_out = VariantFile(vcf_out_filename, 'w', header=vcf_out_hdr, threads=2)

        print("[OUT VCF] will have {} samples from intersecting set".format(set_size))

    delete_previous_outputs(args.out_dir, 'summary.snv' + regionStr, list(samplesDict.keys()))

    # loop over each variant in VCF
    for rec in vcf_in.fetch(rChr, rStart, rEnd):
        #vcf_out.write(rec)

        #if rec.pos < 10684424: continue
        if len(rec.alts) > 1:
            #print("Warning found multiallelic variant")
            continue
        #if ct > 10000:break

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

            # Append subset VFLAGS to INFO field
            rec.info[ "VFLAGS_" + subset ] = ",".join(map(str,vf))

            # Append subset ABHet to INFO field
            rec.info[ "ABHet_" + subset ] = float(abhet)

            # Append VariantType
            rec.info[ "VariantType" ] = vtype

            # count number of vflag(11) for VariantInTargetRatio
            if have_target:
               #vflag_11_ct += (11 not in vf)
               if (11 not in vf):
                  #vflag_11_ct += sum( mi.sa.subsets[subset].values())
                  vflag_11_ct += missing + gt_failed + sum(clean_obs)

        # Append VariantInTargetRatio
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
    print("{0:.2f}".format(end - start))
    write_indiv_summary(prefix_indiv)

    if args.no_output_vcf == False:
        # create index
        time.sleep(1)
        check_output(["tabix", "-f", vcf_out_filename])


def calculate_subgroup_scores(subset, subg, subg_cntl):
    """
    """
    scores = OrderedDict()
    for k in mi.sa.subsets[subset]:
        val = subg[k]
        val_cntl = subg_cntl[k]
        scores['nClean_' + k] = sum(val) # ",".join(map(str,val)),
        scores['Zhet_' + k] = calc_ExcessHet(*val)[0]
        scores['pHWE_' + k] = calc_pHWE(*val_cntl)
        if type(scores['Zhet_' + k]) == float:
            scores['Zhet_' + k] = "{0:.6f}".format(scores['Zhet_' + k])
        if type(scores['pHWE_' + k]) == float:
            scores['pHWE_' + k] = "{0:.6f}".format(scores['pHWE_' + k])

    return scores


def find_s_d(total_obs, samples):
    """
    """
    maf = 0
    if sum(total_obs) > 0:
        maf = (total_obs[1] + 2 * total_obs[2]) / (2 * sum(total_obs))
    if maf <= 0.5:
        # singleton
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
    for k, sm in samples.items():
        if sm['GT'] in {(0,1), (1,0)}:
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


def gather_intersect_fam_vcf_samples(vcf_samples, fam_samples):
    """
    """
    # Method 5 - set, one pass
    for key,sm in vcf_samples.items():
        for subset, sm_set in fam_samples.items():
            if key in sm_set['set']:
                fam_samples[subset]['dict'][key] = sm
    return fam_samples


def check_mendelian_errors(prefix, rec):
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
