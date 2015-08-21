#!/usr/bin/env python
'''
conservation_estimation.py:
Estimate conservation score in miRNAs and targetRNAs.

<Methods>
(1)PhastCons - 46Ways
(2)phyloP - 46way => better

(i)Conservation in TargetRNA seed region
-------- miRNA(8nt_seed)
|||||||| 
-------- TargetRNA(8nt_seed)

(ii)Conservation in candidate target site
---------------------- miRNA
||||||||||||||||||||||
---------------------- TargetRNA


'''
def calc_ave_score(args):
    each_score = list(map(float, args))
    ave_score = int(sum(each_score)/len(each_score)*10000)/10000
    return each_score, ave_score

def calc_region_conservation_ranking(args, cons_score, windows):
    ave_score_list = []
    window_number = len(args)-windows
    for x in range(window_number):
        each_score = list(map(float, args[x:x+windows]))
        ave_score = int((sum(each_score)/len(each_score))*10000)/10000
        ave_score_list.append(ave_score)
    ave_score_list.append(cons_score)
    ave_score_list.sort()
    score_rank = ave_score_list.index(cons_score)/window_number
    return score_rank

def conservation_estimation(mirna_seq, targetrna_seq, mirna_phastcons_dict, mirna_phylop_dict, targetrna_phastcons_dict, targetrna_phylop_dict, tmp_dict):
    seed_length = 8 #windows
    for x in list(tmp_dict.keys()):
        id_infor = x
        data = id_infor.split('||')
        mirna_name_id = data[0]
        targetRNA_id = data[1]
        targetRNA_st = int(data[2]) - 1
        targetRNA_seed_ed = int(data[3])
        mirna_length = len(mirna_seq) #windows
        targetrna_length = len(targetrna_seq)
        targetRNA_ed = targetRNA_st + mirna_length
        
        #Score in 3'UTR region
        targetrna_phastcons_score_all, ave_targetrna_phastcons_score_all = calc_ave_score(targetrna_phastcons_dict[targetRNA_id])
        targetrna_phylop_score_all, ave_targetrna_phylop_score_all = calc_ave_score(targetrna_phylop_dict[targetRNA_id])

        ave_targetrna_phastcons_seed_score = 'NA'
        ave_targetrna_phastcons_score = 'NA'
        ave_targetrna_phylop_seed_score = 'NA'
        ave_targetrna_phylop_score = 'NA'
        diff_phastcons_seed = 'NA'
        diff_phastcons = 'NA'
        diff_phylop_seed = 'NA'
        diff_phylop = 'NA'

        if targetrna_phastcons_dict[targetRNA_id][targetRNA_st:targetRNA_seed_ed] or targetrna_phylop_dict[targetRNA_id][targetRNA_st:targetRNA_seed_ed]:
            #Score in seed matching sites
            targetrna_phastcons_seed_score, ave_targetrna_phastcons_seed_score = calc_ave_score(targetrna_phastcons_dict[targetRNA_id][targetRNA_st:targetRNA_seed_ed])
            targetrna_phastcons_score, ave_targetrna_phastcons_score = calc_ave_score(targetrna_phastcons_dict[targetRNA_id][targetRNA_st:targetRNA_ed])
            targetrna_phylop_seed_score, ave_targetrna_phylop_seed_score = calc_ave_score(targetrna_phylop_dict[targetRNA_id][targetRNA_st:targetRNA_seed_ed])
            targetrna_phylop_score, ave_targetrna_phylop_score = calc_ave_score(targetrna_phylop_dict[targetRNA_id][targetRNA_st:targetRNA_ed])

            #Difference between score in 3'UTR region and seed matching sites
            diff_phastcons_seed = ave_targetrna_phastcons_seed_score - ave_targetrna_phastcons_score_all
            diff_phastcons = ave_targetrna_phastcons_score - ave_targetrna_phastcons_score_all
            diff_phylop_seed = ave_targetrna_phylop_seed_score - ave_targetrna_phylop_score_all
            diff_phylop = ave_targetrna_phylop_score - ave_targetrna_phylop_score_all
        else:
            print('')
            print('WARNINGS: conservation score for ' + targetRNA_id + ' does not exist ...')

        #Score(All_region) - Score(seed_region) - Score(mirna_matching_region) - Diff_Score(seed) - Diff_Score(mirna_matching)
        tmp_dict[x].extend([ave_targetrna_phastcons_score_all, ave_targetrna_phastcons_seed_score, ave_targetrna_phastcons_score, diff_phastcons_seed, diff_phastcons]) #PhastCons_score(Average) => Bad...
        tmp_dict[x].extend([ave_targetrna_phylop_score_all, ave_targetrna_phylop_seed_score, ave_targetrna_phylop_score, diff_phylop_seed, diff_phylop]) #PhyloP(Average) => Good!!

        phastcons_score_rank_seed = 'NA'
        phastcons_score_rank_match = 'NA'
        phylop_score_rank_seed = 'NA'
        phylop_score_rank_match = 'NA'
        
        if targetrna_phastcons_dict[targetRNA_id][targetRNA_st:targetRNA_seed_ed] or targetrna_phylop_dict[targetRNA_id][targetRNA_st:targetRNA_seed_ed]:
            #calculate_conservation_score_ranking_in_seed_region(or miRNA-TargetRNA hybrid region)
            phastcons_score_rank_seed = calc_region_conservation_ranking(targetrna_phastcons_dict[targetRNA_id], ave_targetrna_phastcons_seed_score, seed_length)
            phastcons_score_rank_match = calc_region_conservation_ranking(targetrna_phastcons_dict[targetRNA_id], ave_targetrna_phastcons_score, mirna_length)
            phylop_score_rank_seed = calc_region_conservation_ranking(targetrna_phylop_dict[targetRNA_id], ave_targetrna_phylop_seed_score, seed_length)
            phylop_score_rank_match = calc_region_conservation_ranking(targetrna_phylop_dict[targetRNA_id], ave_targetrna_phylop_score, mirna_length)
        else:
            pass

        tmp_dict[x].extend([phastcons_score_rank_seed, phastcons_score_rank_match, phylop_score_rank_seed, phylop_score_rank_match]) #Larger score is better.

        #print('PhastCons: ',targetrna_phastcons_score)
        #print('PhastCons_ave: ',ave_targetrna_phastcons_score)
        #print('PhyloP: ',targetrna_phylop_score)
        #print('PhyloP_ave: ',ave_targetrna_phylop_score)

    return tmp_dict
