#!/usr/bin/env python
'''
target_site_composition.py:
Calculate target site composition(e.g. AU composition flanking the seed region)

                 -------- miRNA
                 ||||||||
---------------- -------- ---------------- TargetRNA
    upstream       seed      downstream
     (30nt)        (8nt)       (30nt)

<AU content>
(1)seed region
(2)upstream region
(3)downstream region
(4)the ohter 3'UTR region

<AU content difference>
(1)seed region - the ohter 3'UTR region
(2)upstream region - the ohter 3'UTR region
(3)downstream region - the ohter 3'UTR region

[1] Grimson A, Farh KK, Johnston WK, Garrett-Engele P, Lim LP, Bartel DP MicroRNA targeting specificity in mammals: determinants beyond seed pairing. Mol Cell. 2007 Jul 6;27(1):91-105.

'''

def au_contents(args):
    total = 0
    au_number = 0
    for x in args:
        if x == 'A' or x == 'U':
            au_number += 1
        total += 1
    if not total == 0:
        au_content = int(au_number/total*100.0*10000)/10000
    else:
        au_content = 'NA'
    return au_content, au_number

def calc_region_AU_content_ranking(args, cons_score, windows):
    au_score_list = []
    window_number = len(args)-windows
    for x in range(window_number):
        window_seq = args[x:x+windows]
        au_content_window, au_number_window = au_contents(window_seq)
        au_score_list.append(int(au_content_window*10000)/10000)
    au_score_list.append(cons_score)
    au_score_list.sort()
    score_rank = au_score_list.index(cons_score)/window_number
    return score_rank

def target_site_composition(targetrna_seq, tmp_dict, window_size=30):
    seed_length = 8 #windows
    window_length = window_size #windows
    for x in list(tmp_dict.keys()):
        id_infor = x
        data = id_infor.split('||')
        mirna_name_id = data[0]
        targetrna_id = data[1]
        targetrna_st = int(data[2]) - 1
        targetrna_ed = int(data[3])
        seed_up = targetrna_st - window_size
        seed_down = targetrna_ed + window_size

        #region infor
        seed_match = targetrna_seq[targetrna_st:targetrna_ed]
        if seed_up < 0:
            seed_up = 0
        seed_around_up = targetrna_seq[seed_up:targetrna_st] #upstream UTR
        seed_around_down = targetrna_seq[targetrna_ed:seed_down] #downstream UTR
        other_up = targetrna_seq[:seed_up]
        other_down = targetrna_seq[seed_down:]
        total_seq = other_up + other_down

        #AU content in 3'UTR regions
        au_content_utr, au_number_utr = au_contents(total_seq)
        au_content_seed_around_up, au_number_seed_around_up = au_contents(seed_around_up)
        au_content_seed_around_down, au_number_seed_around_down = au_contents(seed_around_down)
        au_content_seed_match, au_number_seed_match = au_contents(seed_match)

        #diff AU content(seed_around vs the other 3'UTR region)
        if au_content_utr == 'NA':
            diff_au_content_seed_around_up = 'NA'
            diff_au_content_seed_around_down = 'NA'
            diff_au_content_seed_match = 'NA'
        else:
            if not au_content_seed_around_up == 'NA':
                diff_au_content_seed_around_up = au_content_seed_around_up - au_content_utr
            else:
                diff_au_content_seed_around_up = 'NA'
            if not au_content_seed_around_down == 'NA':
                diff_au_content_seed_around_down = au_content_seed_around_down - au_content_utr
            else:
                diff_au_content_seed_around_down = 'NA'
            if not au_content_seed_match == 'NA':
                diff_au_content_seed_match = au_content_seed_match - au_content_utr
            else:
                diff_au_content_seed_match = 'NA'

        #calculate_AU-content_ranking_in_seed_region(or miRNA-TargetRNA hybrid region)
        #print(x)
        if not au_content_seed_around_up == 'NA':
            AU_content_rank_upstream = calc_region_AU_content_ranking(targetrna_seq, au_content_seed_around_up, window_length)
        else:
            AU_content_rank_upstream = 'NA'
        if not au_content_seed_around_down == 'NA':
            AU_content_rank_downstream = calc_region_AU_content_ranking(targetrna_seq, au_content_seed_around_down, window_length)
        else:
            AU_content_rank_downstream = 'NA'

        #result(3UTR_AU - seed_up_AU - seed_down_AU - seed_AU - diff_up_AU - diff_down_AU - diff_seed_AU)
        tmp_dict[x].extend([au_content_utr, au_content_seed_around_up, au_content_seed_around_down, au_content_seed_match, diff_au_content_seed_around_up, diff_au_content_seed_around_down, diff_au_content_seed_match])
        tmp_dict[x].extend([AU_content_rank_upstream, AU_content_rank_downstream])

    return tmp_dict
