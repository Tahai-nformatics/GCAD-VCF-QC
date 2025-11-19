#!/usr/bin/env python3

from collections import namedtuple, defaultdict, Counter
import csv
import sys

class Sample:
    def __init__(self, details_dict):
        self.details_dict = details_dict
        self.has_father = False
        self.has_mother = False

    def set_has_father(self):
        self.has_father = True

    def set_has_mother(self):
        self.has_mother = True

    def get_subset(self):
        return self.details_dict.Subset

    def get_subgroup(self):
        return self.details_dict.Subgroup

    def is_control(self):
        """
        AFF: 0 = unknown; 1 = unaffected (controls); 2 = affected (cases)
        """
        return (self.details_dict.AFF == '1')
    
class Sample_multiallelic:
    def __init__(self, details_dict):
        self.details_dict = details_dict
        self.has_father = False
        self.has_mother = False

    def set_has_father(self):
        self.has_father = True

    def set_has_mother(self):
        self.has_mother = True

    def get_subset(self):
        return self.details_dict.Subset

    def get_subgroup(self):
        return self.details_dict.Subgroup

    def is_control(self):
        """
        AFF: 0 = unknown; 1 = unaffected (controls); 2 = affected (cases)
        """
        return (self.details_dict.AFF == '1')
class SampleAnnotation:
    """
    """
    def __init__(self):
        self.sa_collection = dict()
        self.sa_collection_multiallelic = dict()
        self.id_list = set()
        self.subject_list = set()
        #self.good_gt = set()
        self.mi_kids = list()
        self.singletons = list()
        self.private_dbltons = list()
        self.dbltons = list()
        self.subsets = defaultdict(Counter)
        self.subgroups = set()
        self.target_files = set()
        self.subdivide_by_targets = False

    def add_family_sample(self, sample, sample_multiallelic):
        if sample.details_dict.FA in self.id_list:
            sample.set_has_father()
            sample_multiallelic.set_has_father()

        if sample.details_dict.MO in self.id_list:
            sample.set_has_mother()
            sample_multiallelic.set_has_mother()


        # Create a new Sample_multiallelic instance for sa_collection_multiallelic
        self.sa_collection[sample.details_dict.SampID] = sample  # has namedtuple(details_dict)
        self.sa_collection_multiallelic[sample.details_dict.SampID] = sample_multiallelic
    
        # tally
        self.sa_collection[sample.details_dict.SampID].tallySA = dict.fromkeys([-9, (None,None), (0,0), (0,1), (1,0), (1,1),
                                                                                'mend_pair','vp1','vp2',
                                                                                'ti', 'tv', 'non_missing_indel',
                                                                                'ti_wes', 'tv_wes',
                                                                                'singleton', 'p_dblton', 'doubleton',
                                                                                'failing_obs_homo1','failing_obs_het','failing_obs_homo2',
                                                                                'passing_obs_homo1','passing_obs_het','passing_obs_homo2',
                                                                                'missing'],0)

        self.sa_collection_multiallelic[sample.details_dict.SampID].tallySA = dict.fromkeys([-9, (None,None), (0,0), (0,1),(1,0), (1,1),
                                                                                'mend_pair','vp1','vp2',
                                                                                'ti', 'tv', 'non_missing_indel',
                                                                                'ti_wes', 'tv_wes',
                                                                                'singleton', 'p_dblton', 'doubleton',
                                                                                'failing_obs_homo1','failing_obs_het','failing_obs_homo2',
                                                                                'passing_obs_homo1','passing_obs_het','passing_obs_homo2',
                                                                                'missing'],0)


        # DP store
        self.sa_collection[sample.details_dict.SampID].dp_total = 0
        self.sa_collection_multiallelic[sample.details_dict.SampID].dp_total = 0

        # Store subsets and subgroups
        self.subsets[sample.get_subset()][sample.get_subgroup()] += 1
        self.subgroups.add(sample.get_subgroup())

        return

    def add_id(self,indiv_id, subj_id):
        self.id_list.add(indiv_id)
        self.subject_list.add(subj_id)

    def has_father(self,sid):
        return self.sa_collection[sid].has_father

    def has_mother(self,sid):
        return self.sa_collection[sid].has_mother

    def save_good_kid(self,indiv_id):
        if (self.sa_collection[indiv_id].has_father
            or
            self.sa_collection[indiv_id].has_mother):
              self.mi_kids.append(indiv_id)

    def save_good_kid_multiallelic(self,indiv_id):
        if (self.sa_collection_multiallelic[indiv_id].has_father
            or
            self.sa_collection_multiallelic[indiv_id].has_mother):
              self.mi_kids.append(indiv_id)

    def get_mi_kids(self):
        return self.mi_kids

    def get_fam_id(self,indiv_id):
        return self.sa_collection[indiv_id].details_dict.FID

    def get_father(self,sid):
        return self.sa_collection[sid].details_dict.FA

    def get_mother(self,sid):
        return self.sa_collection[sid].details_dict.MO

    def clear_mpairs(self):
        self.mi_kids.clear()

    def tally(self, indiv_id, vsm, failed):
        if failed == 1:
            self.sa_collection[indiv_id].tallySA[ -9 ] += 1
            if vsm['DP'] != None:
                 self.add_dp(indiv_id, vsm['DP'])
            return
        elif failed == -1:
            #if b_or_m == "biallelic":
            # good (passing) genotypes
            self.save_good_kid(indiv_id)
            self.add_dp(indiv_id, vsm['DP'])


        # tally missing and good genotypes
        self.sa_collection[indiv_id].tallySA[ vsm['GT'] ] += 1
    


    def tally_multiallelic(self, indiv_id, vsm, failed):
        if failed == 1:
            self.sa_collection_multiallelic[indiv_id].tallySA[ -9 ] += 1
            if vsm['DP'] != None:
                 self.add_dp(indiv_id, vsm['DP'])
            return
        elif failed == -1:
            #if b_or_m == "biallelic":
            # good (passing) genotypes
            self.save_good_kid_multiallelic(indiv_id)
            self.add_dp_multiallelic(indiv_id, vsm['DP'])

        # tally missing and good genotypes
        if vsm['GT'][0] == vsm['GT'][1]: #Homozygous
            if vsm['GT'] == (0,0): #Ref
                self.sa_collection_multiallelic[indiv_id].tallySA['passing_obs_homo1'] += 1
            else: #Alt
                self.sa_collection_multiallelic[indiv_id].tallySA['passing_obs_homo2'] += 1
        else: #heterozygous
            if 0 in vsm['GT']: #ref het
                self.sa_collection_multiallelic[indiv_id].tallySA['passing_obs_het'] += 1
            else: #alt het
                self.sa_collection_multiallelic[indiv_id].tallySA['passing_obs_homo2'] += 1


    def tally_chrx(self, indiv_id, vsm, failed, chrx_is_multiallelic):
        if failed == 1:
            self.sa_collection_multiallelic[indiv_id].tallySA[ -9 ] += 1
            if vsm['DP'] != None:
                 self.add_dp(indiv_id, vsm['DP'])
            return
        elif failed == -1:
            #if b_or_m == "biallelic":
            # good (passing) genotypes
            self.save_good_kid_multiallelic(indiv_id)
            self.add_dp_multiallelic(indiv_id, vsm['DP'])

        sa_col = self.sa_collection_multiallelic if chrx_is_multiallelic else self.sa_collection
        sex = sa_col[indiv_id].details_dict.SEX

        gt = vsm['GT']

        if gt[0] == gt[1]:  # Homozygous
            key = 'passing_obs_homo1' if gt == (0, 0) else 'passing_obs_homo2'
            sa_col[indiv_id].tallySA[key] += 1
        else:  # Heterozygous
            if 0 in gt:  # Het with reference allele
                key = 'failing_obs_het' if sex == "0" else 'passing_obs_het'
                if sex == "0":
                    sa_col[indiv_id].tallySA['passing_obs_het']
            else:
                key = 'failing_obs_homo2' if sex == "0" else 'passing_obs_homo2'
            sa_col[indiv_id].tallySA[key] += 1


    def tallyIndel(self, ref, alt): #non_missing_indel
        if ref==alt: return
        if len(alt) > 1 or len(ref) > 1:
            self.sa_collection[indiv_id].tallySA['non_missing_indel'] += 1
            return

    def tallyTiTv(self, indiv_id, ref, alt, wes_flag):
        if ref==alt: return
        if len(alt) > 1 or len(ref) > 1:
            return
        
        if ref in {'A', 'G'}:
            if alt in {'A', 'G'}:
                self.sa_collection[indiv_id].tallySA['ti'] += 1
                if wes_flag:
                   self.sa_collection[indiv_id].tallySA['ti_wes'] += 1
            else:
                self.sa_collection[indiv_id].tallySA['tv'] += 1
                if wes_flag:
                   self.sa_collection[indiv_id].tallySA['tv_wes'] += 1
        elif ref in {'C', 'T'}:
            if alt in {'C', 'T'}:
                self.sa_collection[indiv_id].tallySA['ti'] += 1
                if wes_flag:
                   self.sa_collection[indiv_id].tallySA['ti_wes'] += 1
            else:
                self.sa_collection[indiv_id].tallySA['tv'] += 1
                if wes_flag:
                   self.sa_collection[indiv_id].tallySA['tv_wes'] += 1
        
    def tallyTiTv_multiallelic(self, indiv_id, ref, alt, genotype, vtype):

        if vtype == 'MULTI_INDEL':
            return

        for allele_idx in set(genotype):  # Using set to handle homozygous genotypes
            if allele_idx == 0:
                continue
            
            # Get the current alternate allele
            current_alt = alt[allele_idx - 1]  # Note: allele_idx-1 maps to the alt allele
            
            if current_alt == "*": 
                continue

            # If the length of ref and alt are the same, loop through positions
            if len(ref) == len(current_alt):  
                for idx, (a, b) in enumerate(zip(ref, current_alt)):  # Compare corresponding bases
                    if a == b:
                        continue  # No mismatch, skip

                    # If mismatch occurs, check for Ti or Tv
                    if a in {'A', 'G'}:
                        if b in {'A', 'G'}:
                            self.sa_collection_multiallelic[indiv_id].tallySA['ti'] += 1
                        else:
                            self.sa_collection_multiallelic[indiv_id].tallySA['tv'] += 1
                    elif a in {'C', 'T'}:
                        if b in {'C', 'T'}:
                            self.sa_collection_multiallelic[indiv_id].tallySA['ti'] += 1
                        else:
                            self.sa_collection_multiallelic[indiv_id].tallySA['tv'] += 1

            # If the length of ref and alt is different, skip the INDEL cases
            else:
                # For INDELs or other cases where lengths do not match
                continue


    def add_dp(self,indiv_id, dp):
        self.sa_collection[indiv_id].dp_total += dp
    def add_dp_multiallelic(self,indiv_id, dp):
        self.sa_collection_multiallelic[indiv_id].dp_total += dp
    def add_singleton(self, indiv_id):
        self.singletons.append(indiv_id)

    def add_private_dbltons(self, indiv_id):
        self.private_dbltons.append(indiv_id)

    def add_dbltons(self, indiv_id):
        self.dbltons.append(indiv_id)

    def add_target_file(self, trgt_path, subset, target):
        self.target_files.add(trgt_path + ':' + subset + '-' + target)
        self.subdivide_by_targets = True

    def get_targets(self):
        return self.target_files

    def get_divide(self):
        return self.subdivide_by_targets
def createSampleAnnotation(fam):
    """
    """
    delimiter = '\t'
    # Check number of columns
    with open(fam, 'r') as fam_file:
      first_line = fam_file.readline()
    ncol = first_line.count(delimiter) + 1
    if ncol == 1:
        delimiter = ','
        ncol = first_line.count(delimiter) + 1

    if ncol==15:
       SampleFamDetail = namedtuple('SampleFamDetail',['SampID','FID','SubjID','FA','MO','SEX','AFF','AD','AGE','ADSPWGS','Subset','Subgroup','Race_Ethnicity','SeqCtr','ExcludeFromZHet'])

    elif ncol==17:
       SampleFamDetail = namedtuple('SampleFamDetail',['SampID','FID','SubjID','FA','MO','SEX','AFF','AD','AGE','ADSPWGS','Subset','Subgroup','Race_Ethnicity','SeqCtr','ExcludeFromZHet','TargetFile', 'TargetFilePath'])

    else:
      raise TypeError("Unexpected number of columns: %s" % ncol)

    # Read FAM file once to get all sample names
    with open(fam, 'r') as fam_file:
        for sm in map(SampleFamDetail._make, csv.reader(fam_file, delimiter=delimiter)):
            # check values
            sys.tracebacklimit = None
            for s in sm:
                if s in [None, '']:
                    raise ValueError("FAM file should not contain blank values")

            sys.tracebacklimit = 0
            sa.add_id(sm.SampID, sm.SubjID)
            if getattr(SampleFamDetail,'TargetFilePath', None) is not None:
               sa.add_target_file(sm.TargetFilePath, sm.Subset, sm.TargetFile)

    # Re-read FAM file to add in family links
    with open(fam, 'r') as fam_file:
        for sm in map(SampleFamDetail._make, csv.reader(fam_file, delimiter=delimiter)):
            sample = Sample(sm)
            sample_multiallelic = Sample_multiallelic(sm)
            sa.add_family_sample( Sample(sm), sample_multiallelic )
            #print(Sample(sm))
    # remove samples from within subsets having fewer than 5 individuals
    to_delete = list()
    for subset, s_count in sa.subsets.items():
        for k, v in s_count.items():
          if v < 5:
            to_delete.append([subset, k])

    for blk in to_delete:
        print("[FAM] Excluding from HWE {}:{}, too few samples ({})".format(blk[0], blk[1], sa.subsets[ blk[0] ][ blk[1] ] ) )
        #del sa.subsets[ blk[0] ][ blk[1] ]
    return

# Globals
sa = SampleAnnotation()

