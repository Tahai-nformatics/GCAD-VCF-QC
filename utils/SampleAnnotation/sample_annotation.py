#!/usr/bin/env python3

from collections import namedtuple, Counter
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
        return details_dict.Subset

    def get_subgroup(self):
        return details_dict.Subgroup

class SampleAnnotation:
    """
    """
    def __init__(self):
        self.sa_collection = dict()
        self.id_list = set()
        #self.good_gt = set()
        self.mi_kids = list()

    def add_family_sample(self, sample):
        if sample.details_dict.FA in self.id_list:
            sample.set_has_father()

        if sample.details_dict.MO in self.id_list:
            sample.set_has_mother()

        self.sa_collection[sample.details_dict.SampID] = sample # has namedtuple(details_dict)

        # tally
        self.sa_collection[sample.details_dict.SampID].tallySA = Counter()

        # DP store
        self.sa_collection[sample.details_dict.SampID].dp_total = 0
        return

    def add_id(self,indiv_id):
        self.id_list.add(indiv_id)

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

    def tallyTiTv(self, indiv_id, ref, alt):
        if ref==alt: return

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

        """
            if ref == 'G':
                if alt == "A":
                    self.sa_collection[indiv_id].tallySA['ti'] += 1
                else:
                    self.sa_collection[indiv_id].tallySA['tv'] += 1
            elif ref ==  "C":
                if alt == "T":
                    self.sa_collection[indiv_id].tallySA['ti'] += 1
                else:
                    self.sa_collection[indiv_id].tallySA['tv'] += 1
            elif ref ==  "A":
                if alt == "G":
                    self.sa_collection[indiv_id].tallySA['ti'] += 1
                else:
                    self.sa_collection[indiv_id].tallySA['tv'] += 1
            elif ref ==  "T":
                if alt == "C":
                    self.sa_collection[indiv_id].tallySA['ti'] += 1
                else:
                    self.sa_collection[indiv_id].tallySA['tv'] += 1
                    """

    def add_dp(self,indiv_id, dp):
        self.sa_collection[indiv_id].dp_total += dp

def createSampleAnnotation(fam):
    """
    """
    SampleFamDetail = namedtuple('SampleFamDetail',['SampID','FID','SubjID','FA','MO','SEX','AFF','AD','AGE','ADSPWGS','Subset','Subgroup','Race_Ethnicity','SeqCtr','ExcludeFromZHet'])

    with open(fam, 'r') as fam_file:
        for sm in map(SampleFamDetail._make, csv.reader(fam_file, delimiter='\t')):
            sa.add_id(sm.SampID)

    with open(fam, 'r') as fam_file:
        for sm in map(SampleFamDetail._make, csv.reader(fam_file, delimiter='\t')):
            sa.add_family_sample( Sample(sm) )


    return

# Globals
sa = SampleAnnotation()
