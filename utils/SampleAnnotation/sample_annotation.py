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
        self.sa_collection[sample.details_dict.SampID].tallySA = dict.fromkeys([-9, (None,None), (0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(0,7),(0,8),(0,9),(0,10),(0,11),(0,12),(0,13),(0,14),(0,15),(0,16),(1,0),(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(1,7),(1,8),(1,9),(1,10),(1,11),(1,12),(1,13),(1,14),(1,15),(1,16),(2,0),(2,1),(2,2),(2,3),(2,4),(2,5),(2,6),(2,7),(2,8),(2,9),(2,10),(2,11),(2,12),(2,13),(2,14),(2,15),(2,16),(3,0),(3,1),(3,2),(3,3),(3,4),(3,5),(3,6),(3,7),(3,8),(3,9),(3,10),(3,11),(3,12),(3,13),(3,14),(3,15),(3,16),(4,0),(4,1),(4,2),(4,3),(4,4),(4,5),(4,6),(4,7),(4,8),(4,9),(4,10),(4,11),(4,12),(4,13),(4,14),(4,15),(4,16),(5,0),(5,1),(5,2),(5,3),(5,4),(5,5),(5,6),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),(5,13),(5,14),(5,15),(5,16),(6,0),(6,1),(6,2),(6,3),(6,4),(6,5),(6,6),(6,7),(6,8),(6,9),(6,10),(6,11),(6,12),(6,13),(6,14),(6,15),(6,16),(7,0),(7,1),(7,2),(7,3),(7,4),(7,5),(7,6),(7,7),(7,8),(7,9),(7,10),(7,11),(7,12),(7,13),(7,14),(7,15),(7,16),(8,0),(8,1),(8,2),(8,3),(8,4),(8,5),(8,6),(8,7),(8,8),(8,9),(8,10),(8,11),(8,12),(8,13),(8,14),(8,15),(8,16),(9,0),(9,1),(9,2),(9,3),(9,4),(9,5),(9,6),(9,7),(9,8),(9,9),(9,10),(9,11),(9,12),(9,13),(9,14),(9,15),(9,16),(10,0),(10,1),(10,2),(10,3),(10,4),(10,5),(10,6),(10,7),(10,8),(10,9),(10,10),(10,11),(10,12),(10,13),(10,14),(10,15),(10,16),(11,0),(11,1),(11,2),(11,3),(11,4),(11,5),(11,6),(11,7),(11,8),(11,9),(11,10),(11,11),(11,12),(11,13),(11,14),(11,15),(11,16),(12,0),(12,1),(12,2),(12,3),(12,4),(12,5),(12,6),(12,7),(12,8),(12,9),(12,10),(12,11),(12,12),(12,13),(12,14),(12,15),(12,16),(13,0),(13,1),(13,2),(13,3),(13,4),(13,5),(13,6),(13,7),(13,8),(13,9),(13,10),(13,11),(13,12),(13,13),(13,14),(13,15),(13,16),(14,0),(14,1),(14,2),(14,3),(14,4),(14,5),(14,6),(14,7),(14,8),(14,9),(14,10),(14,11),(14,12),(14,13),(14,14),(14,15),(14,16),(15,0),(15,1),(15,2),(15,3),(15,4),(15,5),(15,6),(15,7),(15,8),(15,9),(15,10),(15,11),(15,12),(15,13),(15,14),(15,15),(15,16),(16,0),(16,1),(16,2),(16,3),(16,4),(16,5),(16,6),(16,7),(16,8),(16,9),(16,10),(16,11),(16,12),(16,13),(16,14),(16,15),(16,16),
                                                                                'mend_pair','vp1','vp2',
                                                                                'ti', 'tv', 'non_missing_indel',
                                                                                'ti_wes', 'tv_wes',
                                                                                'singleton', 'p_dblton', 'doubleton',
                                                                                'failing_obs_homo1','failing_obs_het','failing_obs_homo2',
                                                                                'passing_obs_homo1','passing_obs_het','passing_obs_homo2',
                                                                                'missing'],0)



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
            if vsm['DP'] != None:
                 self.add_dp(indiv_id, vsm['DP'])
            return
        elif failed == -1 :
            # good (passing) genotypes
            self.save_good_kid(indiv_id)
            self.add_dp(indiv_id, vsm['DP'])

        # tally missing and good genotypes
        self.sa_collection[indiv_id].tallySA[ vsm['GT'] ] += 1

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
    samplle_id = 0
    subjject_id = 0
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
#        print('opening fam file')
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
            sa.add_family_sample( Sample(sm) )
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
    #print('len is ',len(sa.id_list), len(sa.subject_list)) 
    return

# Globals
sa = SampleAnnotation()
