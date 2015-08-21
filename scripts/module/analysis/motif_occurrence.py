#!/usr/bin/env python
'''
motif_occurence.py:
Calculate probability of the miRNA-binding site occurence

<Methods>
(1)Binomial distribution
(2)Exact probability distribution

<Rules>
WC => Watson-Crick match
GU => G:U Wobble
-------------------------------
seed_length = 8
-------------------------------
8mer: p1-p8 WC match
8mer-1A: p1=A, p2-p8 WC match
7mer-m8: p2-p8 WC match
7mer-m1: p1-p7 WC match
7mer-1A: p1=A, p2-p7 WC match
6mer-m7: p2-p7 WC match
6mer-m8: p3-p8 WC match
-------------------------------
Seed_Grouping
-------------------------------
(1)p1-p8 match: 8-mer
(2)p2-p8 match: 8mer-1A, 7mer-m8
(3)p1-p7 match: 7mer-m1
(4)p2-p7 match: 7mer-1A, 6mer-m7
(5)p3-p8 match: 6mer-m8
-------------------------------

<Reference>
[1] Marín RM, Vanícek J. Efficient use of accessibility in microRNA target prediction. Nucleic Acids Res. 2011 Jan;39(1):19-29.

'''

import itertools
from math import factorial, fsum
from collections import Counter
import utils.setting_utils as utils


def permutations_with_replacement(args, n): #n => string length
    p_list = [''.join(x) for x in list(itertools.product(args,repeat=n))]
    return p_list

#permutations_with_replacement('AUGC',2)
#['AA', 'AU', 'AG', 'AC', 'UA', 'UU', 'UG', 'UC', 'GA', 'GU', 'GG', 'GC', 'CA', 'CU', 'CG', 'CC']

def grouper(iterable, n, fillvalue=None):
    '''Collect data into fixed-length chunks or blocks'''
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)

def Markov_Model(targetrna_seq, nucleotide=['A','U','G','C'], markov_order = 1): #1-order
    n_mer = markov_order + 1
    motifs = permutations_with_replacement(nucleotide,n_mer)
    motif_count_dict = {x:0 for x in motifs}
    for x in range(len(targetrna_seq) - markov_order):
        s_mer = targetrna_seq[x:x+n_mer]
        motif_count_dict[s_mer] += 1
    motif_count = []
    for x in motifs:
        motif_count.append(motif_count_dict[x])
    motifs_group = list(grouper(motifs,len(nucleotide)))
    motif_count_group = list(grouper(motif_count,len(nucleotide)))
    sum_first_nucleotide = list(map(sum, motif_count_group)) #[12,13,21,31]
    output = []
    for x in range(len(nucleotide)): #0, 1, 2, 3
        if int(sum_first_nucleotide[x]) != 0:
            output.extend([motif_count_group[x][y]/float(sum_first_nucleotide[x]) for y in range(len(nucleotide))])
        else:
            output.extend([0 for y in range(len(nucleotide))])
    motif_prob_dict = {motifs[x]:output[x] for x in range(len(motifs))}
    return motif_prob_dict

#P = sum(n!/k!(n!-k!)*p^k*(1-p)^(n-k))
def binomial_distribution(k, n, p):
    combination = factorial(n)/(factorial(k)*factorial(n-k))
    possibility = combination*(p**k)*((1-p)**(n-k))
    return possibility

def cumulative_binomial_distribution(k, n, p):
    frequency = 0.
    for x in range(k+1):
        possibility = binomial_distribution(x, n, p)
        frequency += possibility
    survival = 1.0 - frequency
    if survival < 0:
        survival = 0
    return survival

def seed_grouping(seq):
    seed_group_list = []
    for x in seq:
        if x == '8mer':
            seed_group_list.append('p1_p8_match')
        elif x == '8mer-1A' or x == '7mer-m8':
            seed_group_list.append('p2_p8_match')
        elif x == '7mer-m1':
            seed_group_list.append('p1_p7_match')
        elif x == '7mer-1A' or x == '6mer-m7':
            seed_group_list.append('p2_p7_match')
        elif x == '6mer-m8':
            seed_group_list.append('p3_p8_match')
        else:
            print ('ERROR: motif_type is wrong...')
    seed_group_list = utils.rm_duplicate_list(seed_group_list)
    return seed_group_list

def calc_motif_prob_old_version(seq, two_nt_motif_prob_dict):
    two_nt_motifs = []
    for x in range(len(seq[0])-1):
        seq_list = [seq[y][x:x+2] for y in range(len(seq))]
        seq_list = utils.rm_duplicate_list(seq_list)
        two_nt_motifs.append(seq_list)
    motif_prob = 1.
    #print (two_nt_motifs)
    for x in two_nt_motifs:
        sum_two_nt_motif_prob = 0.
        for each_two_motif in x:
            probability = two_nt_motif_prob_dict[each_two_motif]
            sum_two_nt_motif_prob += probability
            #print (each_two_motif)
            #print (probability)
        motif_prob *= sum_two_nt_motif_prob
        #print ('motif_prob')
        #print (motif_prob)
    #print (motif_prob)
    return motif_prob

def calc_motif_prob(seq, two_nt_motif_prob_dict):
    total_motif_prob = 0.
    each_motif_prob_dict = {} #test
    for motif in seq:
        seq_list = [two_nt_motif_prob_dict[motif[y:y+2]] for y in range(len(motif)-1)] #['AU','AA','AC','AG', ...]
        #print (seq_list)
        motif_prob = 1.
        for y in seq_list:
            motif_prob *= y
        #print (motif_prob)
        each_motif_prob_dict[motif] = motif_prob #test
        total_motif_prob += motif_prob
    #print (total_motif_prob)
    return total_motif_prob, each_motif_prob_dict #test

def motif_occurrence(mirna_seq, targetrna_seq, tmp_dict):
    two_nt_motif_prob_dict = Markov_Model(targetrna_seq)
    #print (two_nt_motif_prob_dict)

    motif_type_need = []
    motif_type_dict = {}
    motif_type_dict['p1_p8_match'] = []
    motif_type_dict['p2_p8_match'] = []
    motif_type_dict['p1_p7_match'] = []
    motif_type_dict['p2_p7_match'] = []
    motif_type_dict['p3_p8_match'] = []

    for x in list(tmp_dict.keys()):
        id_infor = x
        targetrna_motif = tmp_dict[x][3] #NEED TO CHECK!!
        motif_type = tmp_dict[x][5] #NEED TO CHECK!!
        if motif_type == '8mer': #p1_p8_match, p2_p8_match, p1_p7_match, p2_p7_match, p3_p8_match
            motif_type_need.append('p1_p8_match')
            motif_type_dict['p1_p8_match'].append(utils.reverse_complement(targetrna_motif[0:8]))
            motif_type_dict['p2_p8_match'].append(utils.reverse_complement(targetrna_motif[1:8]))
            motif_type_dict['p1_p7_match'].append(utils.reverse_complement(targetrna_motif[0:7]))
            motif_type_dict['p2_p7_match'].append(utils.reverse_complement(targetrna_motif[1:7]))
            motif_type_dict['p3_p8_match'].append(utils.reverse_complement(targetrna_motif[2:8]))
        elif motif_type == '8mer-1A': #p2_p8_match, p2_p7_match, p3_p8_match
            motif_type_need.append('p2_p8_match')
            motif_type_dict['p2_p8_match'].append(utils.reverse_complement(targetrna_motif[1:8]))
            motif_type_dict['p2_p7_match'].append(utils.reverse_complement(targetrna_motif[1:7]))
            motif_type_dict['p3_p8_match'].append(utils.reverse_complement(targetrna_motif[2:8]))
        elif motif_type == '7mer-m8': #p2_p8_match, p2_p7_match, p3_p8_match
            motif_type_need.append('p2_p8_match')
            motif_type_dict['p2_p8_match'].append(utils.reverse_complement(targetrna_motif[1:8]))
            motif_type_dict['p2_p7_match'].append(utils.reverse_complement(targetrna_motif[1:7]))
            motif_type_dict['p3_p8_match'].append(utils.reverse_complement(targetrna_motif[2:8]))
        elif motif_type == '7mer-m1': #p1_p7_match, p2_p7_match
            motif_type_need.append('p1_p7_match')
            motif_type_dict['p1_p7_match'].append(utils.reverse_complement(targetrna_motif[0:7]))
            motif_type_dict['p2_p7_match'].append(utils.reverse_complement(targetrna_motif[1:7]))
        elif motif_type == '7mer-1A': #p2_p7_match
            motif_type_need.append('p2_p7_match')
            motif_type_dict['p2_p7_match'].append(utils.reverse_complement(targetrna_motif[1:7]))
        elif motif_type == '6mer-m7': #p2_p7_match
            motif_type_need.append('p2_p7_match')
            motif_type_dict['p2_p7_match'].append(utils.reverse_complement(targetrna_motif[1:7]))
        elif motif_type == '6mer-m8': #p3_p8_match
            motif_type_need.append('p3_p8_match')
            motif_type_dict['p3_p8_match'].append(utils.reverse_complement(targetrna_motif[2:8]))
        else:
            print ('ERROR: motif_type is wrong...')
    
    motif_type_need = utils.rm_duplicate_list(motif_type_need)
    #motif_type_dict['p1_p8_match'] = utils.rm_duplicate_list(motif_type_dict['p1_p8_match'])
    #motif_type_dict['p2_p8_match'] = utils.rm_duplicate_list(motif_type_dict['p2_p8_match'])
    #motif_type_dict['p1_p7_match'] = utils.rm_duplicate_list(motif_type_dict['p1_p7_match'])
    #motif_type_dict['p2_p7_match'] = utils.rm_duplicate_list(motif_type_dict['p2_p7_match'])
    #motif_type_dict['p3_p8_match'] = utils.rm_duplicate_list(motif_type_dict['p3_p8_match'])
    #print (motif_type_need)
    #print (motif_type_dict['p1_p8_match'])
    #print (motif_type_dict['p2_p8_match'])
    #print (motif_type_dict['p1_p7_match'])
    #print (motif_type_dict['p2_p7_match'])
    #print (motif_type_dict['p3_p8_match'])

    motif_prob_dict = {}
    motif_prob_dict_each = {}

    for x in motif_type_need: #each_type: p1_p8_match, p2_p8_match, p1_p7_match, p2_p7_match, p3_p8_match
        test_seed = utils.rm_duplicate_list(motif_type_dict[x])
        #print (test_seed) #['UGCUUGAA', 'UACUUGAA', 'UAUUUGAG', 'UAUUUGGA']

        #All_existed_motifs_calc
        pos_motif_number = len(motif_type_dict[x])
        total_motif_number = len(targetrna_seq) - len(test_seed[0]) + 1
        motif_prob, each_motif_prob_dict = calc_motif_prob(test_seed,two_nt_motif_prob_dict) #test
        #print ('Motif_number: ' + str(pos_motif_number))
        #print ('Total_motif: ' + str(total_motif_number))
        #print ('Motif_probability: ' + str(motif_prob))
        motif_binom = cumulative_binomial_distribution(pos_motif_number, total_motif_number, motif_prob)
        #print ('Motif_binom_prob: ' + str(motif_binom))
        motif_prob_dict[x] = [pos_motif_number, total_motif_number, motif_prob, motif_binom] #p1_p8_match => [existed_motifs, total_motifs, motif_prob, binom_prob(e.g. 0.00022)]

        #Each_existed_motif_calc
        existed_motif_dict = Counter(motif_type_dict[x])
        #print (existed_motif_dict)
        for i in existed_motif_dict.keys(): #motif => number
            pos_motif_number_each = existed_motif_dict[i]
            motif_prob_each = each_motif_prob_dict[i] #motif => probability
            motif_binom_each = cumulative_binomial_distribution(pos_motif_number_each, total_motif_number, motif_prob_each)
            #print (i)
            #print ('motif_number: ' + str(pos_motif_number_each))
            #print ('motif_prob_each: ' + str(motif_prob_each))
            #print (motif_binom_each)
            motif_prob_dict_each[i] = [pos_motif_number_each, total_motif_number, motif_prob_each, motif_binom_each]

    #print (motif_prob_dict)
    #print (motif_prob_dict_each)

    for x in list(tmp_dict.keys()):
        id_infor = x
        targetrna_motif = tmp_dict[x][3] #NEED TO CHECK
        motif_type = tmp_dict[x][5] #NEED TO CHECK
        if motif_type == '8mer': #p1-p8 match
            seed_group = 'p1_p8_match'
            all_existed_motif_result = motif_prob_dict[seed_group]
            each_existed_motif_result = motif_prob_dict_each[utils.reverse_complement(targetrna_motif[0:8])]
            tmp_dict[x].extend(each_existed_motif_result)
            tmp_dict[x].extend(all_existed_motif_result)
        elif motif_type == '8mer-1A': #p2_p8_match
            seed_group = 'p2_p8_match'
            all_existed_motif_result = motif_prob_dict[seed_group]
            each_existed_motif_result = motif_prob_dict_each[utils.reverse_complement(targetrna_motif[1:8])]
            tmp_dict[x].extend(each_existed_motif_result)
            tmp_dict[x].extend(all_existed_motif_result)
        elif motif_type == '7mer-m8': #p2_p8_match
            seed_group = 'p2_p8_match'
            all_existed_motif_result = motif_prob_dict[seed_group]
            each_existed_motif_result = motif_prob_dict_each[utils.reverse_complement(targetrna_motif[1:8])]
            tmp_dict[x].extend(each_existed_motif_result)
            tmp_dict[x].extend(all_existed_motif_result)
        elif motif_type == '7mer-m1': #p1_p7_match
            seed_group = 'p1_p7_match'
            all_existed_motif_result = motif_prob_dict[seed_group]
            each_existed_motif_result = motif_prob_dict_each[utils.reverse_complement(targetrna_motif[0:7])]
            tmp_dict[x].extend(each_existed_motif_result)
            tmp_dict[x].extend(all_existed_motif_result)
        elif motif_type == '7mer-1A': #p2_p7_match
            seed_group = 'p2_p7_match'
            all_existed_motif_result = motif_prob_dict[seed_group]
            each_existed_motif_result = motif_prob_dict_each[utils.reverse_complement(targetrna_motif[1:7])]
            tmp_dict[x].extend(each_existed_motif_result)
            tmp_dict[x].extend(all_existed_motif_result)
        elif motif_type == '6mer-m7': #p2_p7_match
            seed_group = 'p2_p7_match'
            all_existed_motif_result = motif_prob_dict[seed_group]
            each_existed_motif_result = motif_prob_dict_each[utils.reverse_complement(targetrna_motif[1:7])]
            tmp_dict[x].extend(each_existed_motif_result)
            tmp_dict[x].extend(all_existed_motif_result)
        elif motif_type == '6mer-m8': #p3_p8_match
            seed_group = 'p3_p8_match'
            all_existed_motif_result = motif_prob_dict[seed_group]
            each_existed_motif_result = motif_prob_dict_each[utils.reverse_complement(targetrna_motif[2:8])]
            tmp_dict[x].extend(each_existed_motif_result)
            tmp_dict[x].extend(all_existed_motif_result)
        else:
            print ('ERROR: motif_type is wrong...')

    return tmp_dict #[each_existed_motif_result], [all_existed_motif_result]

    '''
    motif_type_all = []
    motif_type_dict = {}
    for x in list(tmp_dict.keys()):
        motif_targetrna_revcomp = utils.reverse_complement(str(tmp_dict[x][1]))
        types = tmp_dict[x][3]
        motif_type_all.append(tmp_dict[x][3]) #8-mer, 8mer-1A, 7mer-m8, 7mer-m1, 
        motif_type_dict[] = 
    print (motif_type_all)
    motif_type = utils.rm_duplicate_list(motif_type_all)
    print (motif_type)
    seed_group = seed_grouping(motif_type)
    print (seed_group)
    if 'p1_p8_match' in seed_group: #8mer ||||||||
        pass
    elif 'p2_p8_match' in seed_group: #8mer-1A, 7mer-m8, 8mer x|||||||
        pass
    elif 'p1_p7_match, ' in seed_group: #7mer-m1, 8mer |||||||x
        pass
    elif 'p2_p7_match' in seed_group: #7mer-1A, 6mer-m7, 7mer-m1, 8mer-1A, 7mer-m8, 8mer x||||||x
        pass
    elif 'p3_p8_match' in seed_group: #6mer-m8, 8mer-1A, 7mer-m8, 8mer xx||||||
        pass
    else:
        print ('ERROR: seed_group is wrong...')
    '''
    '''
        mirna_infor = x #hsa-miR-26a-5p|MIMAT0000082||NM_005900||936||943
        motif_targetrna = str(tmp_dict[x][1])
        motif_pattern = str(tmp_dict[x][2]) #xx||||||
        motif_targetrna_revcomp = ''
        if (motif_pattern[0] == 'x' or motif_pattern[0] == 'A') and motif_pattern[1] == 'x' and motif_pattern[7] != 'x': #xx||||||, Ax|||||| - 6mer - p1
            motif_targetrna_revcomp = utils.reverse_complement(motif_targetrna[2:8]) #targetrna_motif
            print (motif_targetrna_revcomp)
        elif (motif_pattern[0] == 'x' or motif_pattern[0] == 'A') and motif_pattern[1] != 'x' and motif_pattern[7] == 'x': #x||||||x, A||||||x - 6mer - p2
            motif_targetrna_revcomp = utils.reverse_complement(motif_targetrna[1:7]) #targetrna_motif
        elif (motif_pattern[0] == 'x' or motif_pattern[0] == 'A') and motif_pattern[1] != 'x' and motif_pattern[7] != 'x': #x|||||||, A||||||| - 7mer - p3
            motif_targetrna_revcomp = utils.reverse_complement(motif_targetrna[1:8]) #targetrna_motif
        elif motif_pattern[0] != 'x' and motif_pattern[1] != 'x' and motif_pattern[7] == 'x': #|||||||x - 7mer - p4
            motif_targetrna_revcomp = utils.reverse_complement(motif_targetrna[0:7]) #targetrna_motif
        elif motif_pattern[0] != 'x' and motif_pattern[1] != 'x' and motif_pattern[7] != 'x': #|||||||| - 8mer - p5
            motif_targetrna_revcomp = utils.reverse_complement(motif_targetrna[0:8]) #targetrna_motif
        else:
            print ('ERROR: motif_pattern is wrong...')
    '''
