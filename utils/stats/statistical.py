#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import math

def calc_InbreedingCoeff(aa, ab, bb):
    """
    Inbreeding Coefficient - Although the name Inbreeding Coefficient suggests 
    it is a measure of inbreeding, Inbreeding Coefficient measures the excess heterozygosity 
    at a variant site. It can be used as a proxy for poor mapping (sites that have 
    high Inbreeding Coefficients are typically locations in the genome where the mapping 
    is bad and reads that are in the region mismatch the region because they belong elsewhere). 
    At least 10 samples are required (preferably many more) in order for this annotation to be 
    calculated properly.

    We calculate Inbreeding Coefficient as
            # observed heterozygotes
        1 − ————————————————————————
            # expected heterozygotes

    # expected heterozygotes = 2pq * N
    """
    N = aa+ab+bb
    p = ( (2*aa)+ab )/ (2*N)
    q = 1 - p
    return 1 - ( ab / (2*p*q*N) )

def calc_ExcessHet(aa, ab, bb):
    """
    ExcessHet - Phred-scaled p-value for exact test of excess heterozygosity
    This annotation estimates excess heterozygosity in a population of samples. 
    It is related to but distinct from InbreedingCoeff, which estimates evidence 
    for inbreeding in a population. ExcessHet scales more reliably to large cohort sizes.

    p=freq of allele a in population
    q=1-p
    """
    N = aa+ab+bb

    result = '.'
    maf = 0
    if N:
        maf = ((2 * bb) + ab)/(2*N)
        #p =   ((2 * aa) + ab)/(2*N)
        #q = 1 - p

    if maf > 0.5: maf = 1 - maf

    if (ab == N):
        result = '.'
    elif ( (maf >0) and (ab >0) ):
        hetExpct = 2 * maf * ( 1 - maf)
        hetObs   = ab / N
        t  = hetObs - hetExpct
        chisq = (N * t * t) / ( hetObs * (1 - hetObs) )
        result = math.sqrt(chisq)
        if t < 0:
            result = -result

    return [result, maf]

def calc_ExcessHet_OLD(aa, ab, bb):
    from subprocess import check_output
    zhet, hetz_maf = (check_output(
                                ["perl",
                                    "-e",
                                    "use lib '/home/ottov/gcad-qc-pipeline';use ExcessHeterozygosity;@a=ExcessHeterozygosity::excess_het(%d,%d,%d);print join('\t',@a);" %
                                    (aa, ab, bb)])
                            .decode('ascii').split('\t'))
    return [zhet, hetz_maf]

def calc_pHWE(aa,ab,bb):
    """
    """

    if (aa < 0
        or ab < 0
        or bb < 0): return -1

    if(aa < bb):
        obs_homr = aa
        obs_homc = bb
    else:
        obs_homr = bb
        obs_homc = aa


    rare_copies = 2 * obs_homr + ab
    genotypes = obs_homc + obs_homr + ab

    if genotypes == 0: return 0
    elif genotypes < 0: return -1

    mid = int(rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes))

    if ((rare_copies & 1) ^ (mid & 1)):
        mid += 1

    het_probs = [0.0 for x in range(rare_copies + 1)]
    total = het_probs[mid] = 1.0

    # Bottom mid
    curr_hets = mid
    curr_homr = (rare_copies - mid) / 2
    curr_homc = genotypes - curr_hets - curr_homr

    # Bottom scan
    while curr_hets > 1:
        het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
        total += het_probs[curr_hets - 2]
        curr_homr += 1
        curr_homc += 1
        curr_hets -= 2

    # Top mid
    curr_hets = mid;
    curr_homr = (rare_copies - mid) / 2;
    curr_homc = genotypes - curr_hets - curr_homr;

    # Top scan
    while curr_hets <= (rare_copies - 2):
        het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
        total += het_probs[curr_hets + 2]
        curr_homr -= 1
        curr_homc -= 1
        curr_hets += 2

    for i in range(rare_copies + 1):
        het_probs[i] /= total

    p_hwe = 0.0

    for i in range(rare_copies + 1):
        if (het_probs[i] > het_probs[ab]): continue

        p_hwe += het_probs[i]

    if p_hwe > 1: p_hwe = '.'

    return p_hwe

def calc_pHWE_OLD(obs_hom1, obs_hets, obs_hom2):
    """ 
    """
    from subprocess import check_output

    p_hwe = (check_output(["perl",
                            "-e", 
                            "use lib '/home/ottov/gcad-qc-pipeline';use HardyWeinberg;print HardyWeinberg::snphwe(%d,%d,%d);" %
                                (obs_hets, obs_hom1, obs_hom2)]).decode('ascii'))

    if p_hwe != '.':p_hwe = "{:.4f}".format(float(p_hwe))
    if p_hwe == 1: p_hwe='.'

    return p_hwe

def calc_ABHet(ad, dp):
    """
    AlleleBalanceHet - Allele Balance for heterozygous calls (ref/(ref+alt))
    {# REF reads from heterozygous samples}/{# REF + ALT reads from heterozygous samples}

    """
    if dp:
        return ad/dp
    else:
        return 0
