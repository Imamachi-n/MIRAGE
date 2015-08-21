#!/usr/bin/env python

"""
MIRAGE estimate: Comprehensive miRNA target prediction pipeline.

Created by Naoto Imamachi on 2015-04-23.
Copyright (c) 2015 Naoto Imamachi. All rights reserved.
Updated and maintained by Naoto Imamachi since Apr 2015.

"""

import os, sys, time
import yaml
import itertools
import shelve
import utils.setting_utils as utils
from parameter.common_parameters import common_parameters
from parameter.mirna_seed_infor import mirna_seed_infor
from module.analysis.seed_match import find_mirna_target_candidates
from module.analysis.thermo_calc import calc_thermo
from module.analysis.motif_occurrence import motif_occurrence
from module.analysis.conservation_estimation import conservation_estimation
from module.analysis.target_site_composition import target_site_composition

#test
from module.analysis.seed_match_rev import find_mirna_subtarget_candidates

p = utils.Bunch(common_parameters)
seed_infor = utils.Bunch(mirna_seed_infor)
#print (seed_infor.seed_length)

flg_find_mirna_target_candidates = 0

def progress_bar(per, counter, barlen=50):
    perb = int(per/(100.0/barlen))
    s = '\r'
    s += '['
    s += '*' * perb
    s += ' ' * (barlen - perb)
    s += ']'
    s += ' ' + (str(int(per)) + '%' + ' - ' + str(counter) + ' pairs').rjust(4)
    sys.stdout.write(s)

def run_tasks(args, function):
    for x in args:
        result_data = function(result_data)
        if not result_data == None:
            seed_match_db.update(result_data)

def get_sequence(mirna_id, targetrna_id):
    mirna_seq = mirna_dict[mirna_id]
    targetrna_seq = targetrna_dict[targetrna_id]
    return mirna_seq, targetrna_seq

def run_log(comment, step):
    global flg_find_mirna_target_candidates
    if flg_find_mirna_target_candidates == step:
        utils.now_time(comment)
        flg_find_mirna_target_candidates += 1

def run_result(result_dict):
    if not result_dict: #return 'None' if there is not 'tmp_dict' directory data
        return
    else:
        return result_dict

def detect_seed_match(mirna_id, targetrna_id):
    '''
    1_Find_mirna_target_candidates
    Start find_mirna_target_candidates module
    
    '''
    #seq_data
    mirna_seq, targetrna_seq = get_sequence(mirna_id, targetrna_id)
    targetrna_seq_revcomp = utils.reverse_complement(targetrna_seq)
    #print (mirna_seq)
    #print (targetrna_seq)
    #print (targetrna_seq_revcomp)
    #parameters - seed
    mirna_start_pairing = ''
    seed_length = ''
    allowed_gu_wobbles = []
    allowed_mismatches = []
    if hasattr(seed_infor,'MIRNA_START_PAIRING'):
        mirna_start_pairing = seed_infor.MIRNA_START_PAIRING
    else:
        print ('ERROR: MIRNA_START_PAIRING parameter does not exist in module.analysis.mirna_seed.py file')
        sys.exit(1)
    if hasattr(seed_infor,'SEED_LENGTH'):
        seed_length = seed_infor.SEED_LENGTH
    else:
        print ('ERROR: SEED_LENGTH parameters do not exist in module.analysis.mirna_seed.py file')
        sys.exit(1)
    if hasattr(seed_infor,'ALLOWED_GU_WOBBLES'):
        allowed_gu_wobbles = seed_infor.ALLOWED_GU_WOBBLES
    else:
        print ('ERROR: ALLOWED_GU_WOBBLES parameters do not exist in module.analysis.mirna_seed.py file')
        sys.exit(1)
    if hasattr(seed_infor,'ALLOWED_MISMATCHES'):
        allowed_mismatches = seed_infor.ALLOWED_MISMATCHES
    else:
        print ('ERROR: ALLOWED_MISMATCHES parameters do not exist in module.analysis.mirna_seed.py file')
        sys.exit(1)

    tmp_dict = {}
    #run_log("Finding seed matches and calculating motif density in targetRNA sequences...", 0)
    tmp_dict = find_mirna_target_candidates(mirna_id,mirna_seq,targetrna_id,targetrna_seq_revcomp,mirna_start_pairing,seed_length,allowed_gu_wobbles,allowed_mismatches) # => list()
    return run_result(tmp_dict)

def calc_thermo_func(mirna_id, targetrna_id, tmp_dict):
    '''
    2_Calculate_thermodynamic_stability
    Start find_mirna_target_candidates module

    '''
    #seq_data
    mirna_seq, targetrna_seq = get_sequence(mirna_id, targetrna_id)

    #run_log("Calculating thermodynamic stability between miRNA and targetRNA sequences...", 1)
    targetrna_range = 70
    result_dict = calc_thermo(mirna_seq, targetrna_seq, targetrna_range, tmp_dict)
    return run_result(result_dict)

def motif_occurrence_func(mirna_id, targetrna_id, tmp_dict):
    '''
    3_Probalility_of_miRNA_binding_site_occurrence

    '''
    #seq_data
    mirna_seq, targetrna_seq = get_sequence(mirna_id, targetrna_id)

    #run_log("Calculating miRNA-binding sites occurence in targetRNAs...", 2)
    result_dict = motif_occurrence(mirna_seq, targetrna_seq, tmp_dict)
    return run_result(result_dict)

def conservation_estimation_func(mirna_id, targetrna_id, tmp_dict):
    '''
    4_Estimate_conservation_of_miRNAs_and_targetRNAs

    '''
    #seq_data
    mirna_seq, targetrna_seq = get_sequence(mirna_id, targetrna_id)

    #run_log("Estimating conservation score in miRNAs and targetRNAs...", 3)
    result_dict = conservation_estimation(mirna_seq, targetrna_seq, mirna_phastcons_dict, mirna_phylop_dict, targetrna_phastcons_dict, targetrna_phylop_dict, tmp_dict)
    return run_result(result_dict)

def target_site_composition_func(mirna_id, targetrna_id, tmp_dict):
    '''
    5_Calculate_target_site_composition

    '''
    #seq_data
    mirna_seq, targetrna_seq = get_sequence(mirna_id, targetrna_id)

    #run_log("Calculating target site composition...", 4)
    result_dict = target_site_composition(targetrna_seq, tmp_dict)
    return run_result(result_dict)
    
def detect_rev_seed_match(mirna_id, targetrna_id):
    '''
    X1_seed_match_rev

    '''
    mirna_seq, targetrna_seq = get_sequence(mirna_id, targetrna_id)
    targetrna_seq_revcomp = utils.reverse_complement(targetrna_seq)
    tmp_dict = find_mirna_subtarget_candidates(mirna_id,mirna_seq,targetrna_id,targetrna_seq_revcomp) # => list()
    return run_result(tmp_dict)

###MAIN###
utils.now_time("MIRAGE estimate is starting...")
mirna_dict = utils.load_fasta(p.MIRNA_FASTA_PATH)
targetrna_dict = utils.load_fasta(p.TARGETRNA_FASTA_PATH)

'''#shelve
#shelve_file
###Save_file
shelve_path = utils.get_absolute_path('./seed_match.db')
if os.path.isfile(shelve_path): #if shelve_file exists, it'll be removed.
    os.remove(shelve_path)
seed_match_db = shelve.open('./seed_match.db')
'''#shelve

###Conservation_files
mirna_phastcons_path = utils.get_absolute_path('../data/PhastCons46Ways/phastCons46way_miRBase_v21_hg38Tohg19_for_MIRAGE.db')
mirna_phylop_path = utils.get_absolute_path('../data/PhyloP/phyloP46way_miRBase_v21_hg38Tohg19_for_MIRAGE.db')
#targetrna_phastcons_path = utils.get_absolute_path('../data/PhastCons46Ways/phastCons46way_Refseq_for_MIRAGE.db')
#targetrna_phylop_path = utils.get_absolute_path('../data/PhyloP/phyloP46way_Refseq_for_MIRAGE.db')
targetrna_phastcons_path = utils.get_absolute_path('../data/PhastCons46Ways/phastCons46way_Refseq_for_MIRAGE_CDS.db')
targetrna_phylop_path = utils.get_absolute_path('../data/PhyloP/phyloP46way_Refseq_for_MIRAGE_CDS.db')
mirna_phastcons_dict = shelve.open(mirna_phastcons_path)
mirna_phylop_dict = shelve.open(mirna_phylop_path)
targetrna_phastcons_dict = shelve.open(targetrna_phastcons_path)
targetrna_phylop_dict = shelve.open(targetrna_phylop_path)

#mirna_targetrna_list = list(itertools.product(['hsa-miR-26a-5p|MIMAT0000082','hsa-miR-143-3p|MIMAT0000435'],['NM_005900','NM_002658']))
#mirna_targetrna_list = list(itertools.product(['hsa-miR-155-5p|MIMAT0000646'],['NM_005375']))
#mirna_targetrna_list = list(itertools.product(['hsa-miR-155-5p|MIMAT0000646'],['NM_005900']))
#mirna_targetrna_list = list(itertools.product(['hsa-miR-26a-5p|MIMAT0000082'],['NM_005900']))
#mirna_targetrna_list = list(itertools.product(['hsa-miR-155-5p|MIMAT0000646'],['NM_004781']))

#mirna_targetrna_list = list(itertools.product(mirna_dict.keys(),targetrna_dict.keys()))
#mirna_targetrna_list = list(itertools.product(['hsa-miR-155-5p|MIMAT0000646'],targetrna_dict.keys()))
mirna_targetrna_list = list(itertools.product(['hsa-miR-124-3p|MIMAT0000422'],targetrna_dict.keys()))
#mirna_targetrna_list = list(itertools.product(['hsa-miR-155-5p|MIMAT0000646'],['NM_001288811']))


#'''
#test
test_output = open('./miRNA_TargetRNA_list.txt','w')

for line in mirna_targetrna_list:
    print("\t".join(line),end="\n",file=test_output)

test_output.close()
#'''

'''#rev_seed
#test - rev_seed_matching
rev_seed_output = open('./MIRAGE_result_rev_seed_matching.txt','w')
#header
print('miRNA_name_id','TargetRNA_id','TargetRNA_seed_st','TargetRNA_seed_ed','TargetRNA_length',
      'miRNA_match','TargetRNA_match','miRNA_seed_match','TargetRNA_seed_match',
      'seed_match', sep="\t",end="\n",file=rev_seed_output)

mirna_targetrna_number = len(mirna_targetrna_list)
utils.now_time("Testing miRNA-TargetRNA pairs: " + str(mirna_targetrna_number) + ' pairs')

counter = 0
for x in mirna_targetrna_list:
    mirna_id, targetrna_id = x
    #print(x[1])
    counter += 1
    result_data1, key_list = detect_rev_seed_match(mirna_id, targetrna_id)
    #print(result_data1,key_list)
    if result_data1 == 0:
        continue
    else:
        #result_data2 = calc_thermo_func(mirna_id, targetrna_id, result_data1)
        for key in key_list:
            mirna_infor = key
            mirna_data = mirna_infor.split('||')
            mirna_name = mirna_data[0]
            targetrna_name = mirna_data[1]
            targetrna_length = str(len(targetrna_dict[targetrna_name]))
            mirna_st = mirna_data[2]
            mirna_ed = mirna_data[3]
            targetrna_infor = map(str,result_data1[key])
            targetrna_output = "\t".join(targetrna_infor)
            #targetrna_output_flatten = utils.flatten(targetrna_output)
            mirna_output = "\t".join([mirna_name, targetrna_name, mirna_st, mirna_ed, targetrna_length])
            print (mirna_output, targetrna_output, sep="\t", end="\n", file=rev_seed_output)
    par = counter/mirna_targetrna_number*100
    progress_bar(par, counter)
progress_bar(100, counter)
print('')

rev_seed_output.close()
'''#rev_seed

#result_output
result_output = open('./mirage_output.result','w')

#header
print('miRNA_name_id','TargetRNA_id','TargetRNA_seed_st','TargetRNA_seed_ed','TargetRNA_length',
      'miRNA_match','TargetRNA_match','miRNA_seed_match','TargetRNA_seed_match',
      'seed_match','seed_type','GU_wobbles', #'subseed_1','subseed_2',
      'site_density_simple','site_density_MuTaMe',
      'str_mfe_seed', 'mfe_seed', 'str_ens_seed', 'ens_seed', 'delta_G_binding_seed',
      'str_mfe', 'mfe', 'str_ens', 'ens', 'delta_G_binding',
      '3_pairing_contrib_mfe', '3_pairing_contrib_ens', '3_pairing_contrib_delta_G_binding',
      'seed_match_motif_number_each', 'total_motif_number_each', 'motif_prob_each', 'motif_binom_prob_each',
      'seed_match_motif_number_all', 'total_motif_number_all', 'motif_prob_all', 'motif_binom_prob_all',
      'PhastCons_score_ave_all', 'PhastCons_score_ave_seed', 'PhastCons_score_ave_match', 'diff_PhastCons_seed', 'diff_PhastCons_match',
      'PhyloP_score_ave_all', 'PhyloP_score_ave_seed', 'PhyloP_score_ave_match', 'diff_PhyloP_seed', 'diff_PhyloP_match',
      'PhastCons_score_rank_seed', 'PhastCons_score_rank_match', 'PhyloP_score_rank_seed', 'PhyloP_score_rank_match',
      'AU-content_utr', 'AU-content_seed_around_up', 'AU-content_seed_around_down', 'AU-content_seed_match',
      'diff_AU-content_seed_around_up', 'diff_AU-content_seed_around_down', 'diff_AU-content_seed_match',
      'AU_content_rank_upstream', 'AU_content_rank_downstream'
       ,sep="\t", end="\n", file=result_output)


mirna_targetrna_number = len(mirna_targetrna_list)
utils.now_time("Testing miRNA-TargetRNA pairs: " + str(mirna_targetrna_number) + ' pairs')
progress_bar(0, 0)
counter = 0
for x in mirna_targetrna_list:
    mirna_id, targetrna_id = x
    #print(x[1])
    counter += 1
    result_data1, key_list = detect_seed_match(mirna_id, targetrna_id)
    if result_data1 == 0:
        continue
    else:
        #result_data2 = calc_thermo_func(mirna_id, targetrna_id, result_data1)
        #result_data3 = motif_occurrence_func(mirna_id, targetrna_id, result_data2)
        #result_data4 = conservation_estimation_func(mirna_id, targetrna_id, result_data3)
        #result_data5 = target_site_composition_func(mirna_id, targetrna_id, result_data4)
        #seed_match_db.update(result_data5)
        for key in key_list:
            mirna_infor = key
            mirna_data = mirna_infor.split('||')
            mirna_name = mirna_data[0]
            targetrna_name = mirna_data[1]
            targetrna_length = str(len(targetrna_dict[targetrna_name]))
            mirna_st = mirna_data[2]
            mirna_ed = mirna_data[3]
            targetrna_infor = map(str,result_data1[key]) #mod
            targetrna_output = "\t".join(targetrna_infor)
            #targetrna_output_flatten = utils.flatten(targetrna_output)
            mirna_output = "\t".join([mirna_name, targetrna_name, mirna_st, mirna_ed, targetrna_length])
            print (mirna_output, targetrna_output, sep="\t", end="\n", file=result_output)
    par = counter/mirna_targetrna_number*100
    progress_bar(par, counter)
progress_bar(100, counter)
print('')

'''#shelve
for x in list(seed_match_db.keys()):
    mirna_infor = x
    mirna_data = mirna_infor.split('||')
    mirna_name = mirna_data[0]
    targetrna_name = mirna_data[1]
    targetrna_length = str(len(targetrna_dict[targetrna_name]))
    mirna_st = mirna_data[2]
    mirna_ed = mirna_data[3]
    targetrna_infor = map(str,seed_match_db[x])
    targetrna_output = "\t".join(targetrna_infor)
    #targetrna_output_flatten = utils.flatten(targetrna_output)
    mirna_output = "\t".join([mirna_name, targetrna_name, mirna_st, mirna_ed, targetrna_length])
    print (mirna_output, targetrna_output, sep="\t", end="\n", file=result_output)

result_output.close()
'''#shelve

#Finished
utils.now_time("MIRAGE estimate was successfully finished!!")
