#!/usr/bin/env python3

from collections import namedtuple, defaultdict, Counter
import csv

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

class SampleAnnotation:
    """
    """
    def __init__(self):
        self.sa_collection = dict()
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

    def add_family_sample(self, sample):
        if sample.details_dict.FA in self.id_list:
            sample.set_has_father()

        if sample.details_dict.MO in self.id_list:
            sample.set_has_mother()

        self.sa_collection[sample.details_dict.SampID] = sample # has namedtuple(details_dict)

        # tally
        self.sa_collection[sample.details_dict.SampID].tallySA = dict.fromkeys([-9, (None,None), (0,0), (0,1), (1,0), (1,1),
                                                                                'mend_pair','vp1','vp2',
                                                                                'ti', 'tv', 'non_missing_indel',
                                                                                'singleton', 'p_dblton', 'doubleton'],0)
        # DP store
        self.sa_collection[sample.details_dict.SampID].dp_total = 0

        # store the subsets-subgroups
#        self.subsets[ sample.get_subset() ].add(sample.get_subgroup())
        self.subsets[ sample.get_subset() ][ sample.get_subgroup() ] += 1

        # store the subgroups
        self.subgroups.add( sample.get_subgroup() )

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
            return
        elif failed == -1 :
            # good (passing) genotypes
            self.save_good_kid(indiv_id)
            self.add_dp(indiv_id, vsm['DP'])

        # tally missing and good genotypes
        self.sa_collection[indiv_id].tallySA[ vsm['GT'] ] += 1

    def tallyTiTv(self, indiv_id, ref, alt):  #non_missing_indel
        if ref==alt: return
        if len(alt) > 1 or len(ref) > 1:
            self.sa_collection[indiv_id].tallySA['non_missing_indel'] += 1
            return

        if ref in {'A', 'G'}:
            if alt in {'A', 'G'}:
                self.sa_collection[indiv_id].tallySA['ti'] += 1
            else:
                self.sa_collection[indiv_id].tallySA['tv'] += 1
        elif ref in {'C', 'T'}:
            if alt in {'C', 'T'}:
                self.sa_collection[indiv_id].tallySA['ti'] += 1
            else:
                self.sa_collection[indiv_id].tallySA['tv'] += 1

    def add_dp(self,indiv_id, dp):
        self.sa_collection[indiv_id].dp_total += dp

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

    # Check number of columns
    with open(fam, 'r') as fam_file:
      first_line = fam_file.readline()

    ncol = first_line.count('\t') + 1

    if ncol==15:
       SampleFamDetail = namedtuple('SampleFamDetail',['SampID','FID','SubjID','FA','MO','SEX','AFF','AD','AGE','ADSPWGS','Subset','Subgroup','Race_Ethnicity','SeqCtr','ExcludeFromZHet'])

    elif ncol==17:
       SampleFamDetail = namedtuple('SampleFamDetail',['SampID','FID','SubjID','FA','MO','SEX','AFF','AD','AGE','ADSPWGS','Subset','Subgroup','Race_Ethnicity','SeqCtr','ExcludeFromZHet','TargetFile', 'TargetFilePath'])

    else:
      raise TypeError("Unexpected number of columns: %s" % ncol)

    # Read FAM file once to get all sample names
    with open(fam, 'r') as fam_file:
        for sm in map(SampleFamDetail._make, csv.reader(fam_file, delimiter='\t')):
            sa.add_id(sm.SampID, sm.SubjID)
            if sm.TargetFilePath:
               sa.add_target_file(sm.TargetFilePath, sm.Subset, sm.TargetFile)

    # Re-read FAM file to add in family links
    with open(fam, 'r') as fam_file:
        for sm in map(SampleFamDetail._make, csv.reader(fam_file, delimiter='\t')):
            sa.add_family_sample( Sample(sm) )

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
