#!/usr/bin/env python
'''
Seed_match.py:
Check seed matching between miRNA and TargetRNA

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
-------------------------------
seed_length > 8
-------------------------------
9mer: p1-p9 WC match
9mer-1A: p1=A, p2-p9 WC match
10mer: p1-p10 WC match
10mer-1A: p1=A, p2-p10 WC match
...
-------------------------------
Default:
8mers => allowed 2 G:U Wobbles
7mers => allowed 1 G:U Wobble
6mers => allowed 0 G:U Wobble

<Exception>
CLASH or PAR-CLIP analysis based on NGS revealed that there are non-canonical miRNA-targetRNA interactions[1].
-Frequent symmetrical loop sequences
G  U  U  A => miRNA
:  :  x  x
U  G  U  G => TargetRNA

-Bulged pivot
Target mRNAs sometimes have a pivot nucleotide.

<Reference>
[1] Hausser J, Zavolan M. Identification and consequences of miRNA-target interactions--beyond repression of gene expression. Nat Rev Genet. 2014 Sep;15(9):599-612.

'''

from module.analysis.site_density import site_density_simple, site_density_MuTaMe

def is_gu_wobbles(mirna_nuc,targetrna_nuc):
    #check if mirna-targetrna pair is GU-wobble or not.
    if (mirna_nuc == 'G' and targetrna_nuc == 'A') or (mirna_nuc == 'U' and targetrna_nuc == 'C'):
        return True
    else:
        return False

def find_seed_pairings(mirna_seed,targetrna_seed,seed_length=8, allowed_gu_wobbles={6:0,7:1,8:2}, allowed_mismatches={6:0,7:0,8:0}):
    #base pairing - 1bp resolution
    seed_pairs_list = [[mirna_seed[x],targetrna_seed[x]] for x in range(0,seed_length)] #[[miRNA_seed1,TargetRNA_seed1]...]

    #param
    s_mer = ''
    s_class = ''

    fail_criteria = seed_length - 5 #Default:3 => seed_seq < 6
    seed_match = ''
    match_count = 0
    mismatch_count = 0
    wobble_count = 0

    #seed_pairing
    for x in seed_pairs_list:
        mirna_nuc = x[0]
        targetrna_nuc = x[1]
        if mirna_nuc == targetrna_nuc: #WC match
            seed_match += '|'
            match_count += 1
        elif is_gu_wobbles(mirna_nuc,targetrna_nuc): #GU-wobble match
            seed_match += ':'
            wobble_count += 1
            if wobble_count >= fail_criteria:
                return 0,0,1
        else: #mismatch
            seed_match += 'x'
            mismatch_count += 1
            if mismatch_count >= fail_criteria:
                return 0,0,2

    #first_nucleotide_in_TargetRNA => A or not
    if (seed_pairs_list[0][1] == 'U' and seed_match[0] == 'x'):
        seed_match = 'A' + seed_match[1:]
        #print(seed_match)

    #print(seed_match)
    #print(seed_pairs_list[0][1])

    #check_followed_wobbles
    if (wobble_count == 2 and seed_match[0] == 'A'): #Not allowed 2 GU-wobbles & A nucleotide in first bp in TargetRNA
        return 0,0,3

    #check_seed_core(3-7bp)
    seed_3_7 = seed_match[2:7]
    if 'x' in seed_3_7:
        return 0,0,4

    #check_seed_side(2bp,8bp)
    seed_1 = seed_match[0]
    seed_2 = seed_match[1]
    seed_8 = seed_match[7]
    if (seed_2 == 'x' and seed_8 == 'x'): #mismatch in seed core
        return 0,0,5
    elif (seed_2 == 'x' and seed_8 != 'x'): #6mer-m8
        if seed_8 == '|' and wobble_count == 0: #6mer-m8 ['|x||||||',':x||||||','Ax||||||','xx||||||'] / allowed 0 GU-wobble
            s_mer = '6mer'
            s_class = '-m8'
        else:
            return 0,0,6
    elif (seed_2 != 'x' and seed_8 == 'x'): #6mer-m7|7mer-1A|7mer-m1
        if seed_1 == 'x': #6mer-m7 ['x||||||x']
            if wobble_count == 0:
                s_mer = '6mer'
                s_class = '-m7'
            else: #6mer => allowed 0 GU-wobble
                return 0,0,7
        elif seed_1 == 'A': #7mer-1A ['A||||||x','A:|||||x','A|:||||x','A||:|||x','A|||:||x','A||||:|x','A|||||:x']
            if wobble_count <= 1:
                s_mer = '7mer'
                s_class = '-1A'
            else:
                return 0,0,8
        elif wobble_count <= 1: #7mer-m1 ['|||||||x',':||||||x','|:|||||x','||:||||x','|||:|||x','||||:||x','|||||:|x','||||||:x']
            s_mer = '7mer'
            s_class = '-m1'
        else:
            return 0,0,9
    elif (seed_2 != 'x' and seed_8 != 'x'): #7mer-8m|8mer-1A|8mer
        if seed_1 == 'x': #7mer-8m ['x|||||||','x:||||||','x|:|||||','x||:||||','x|||:|||','x||||:||','x|||||:|','x||||||:'] / allowed 1 GU-wobble
            if wobble_count <= 1:
                s_mer = '7mer'
                s_class = '-m8'
            else:
                return 0,0,10
        else: #8mer-1A|8mer
            if seed_1 == 'A': #8mer-1A ['A:|||||||','A|:||||||','A||:|||||','A|||:||||','A||||:|||','A|||||:||','A||||||:|','A|||||||:',]
                s_mer = '8mer'
                s_class = '-1A'
            else: #8mer ['||||||||',':|||||||','|:||||||','||:|||||','|||:||||','||||:|||','|||||:||','||||||:|','|||||||:',
                        #'::||||||',':|:|||||',':||:||||',':|||:|||',':||||:||',':|||||:|',':|||||:|',':||||||:',
                        #'|::|||||','|:|:||||','|:||:|||','|:|||:||','|:||||:|','|:|||||:',
                        #'||::||||','||:|:|||','||:||:||','||:|||:|','||:||||:',
                        #'|||::|||','|||:|:||','|||:||:|','|||:|||:'
                        #'||||::||','||||:|:|','||||:||:',
                        #'|||||::|','|||||:|:'
                        #'||||||::']
                s_mer = '8mer'
                s_class = ''
    else:
        return 0,0,11

    s_type = s_mer + s_class
    return seed_match, s_type, wobble_count

    #return seed_match


def find_mirna_target_candidates(mirna_id, mirna_seq, targetrna_id, targetrna_revcomp_seq, mirna_start_pairing=1, seed_length=8, allowed_gu_wobbles={6:0,7:1,8:2}, allowed_mismatches={6:0,7:0,8:0}):
    """Searches for seed(s) in the target sequence.
        Default_param:
            mirna_start_pairing: 1
            seed_length: 8
            allowed_gu_wobbles: {6:0,7:1,8:2}
            allowed_mismatches: {6:0,7:0,8:0}
    """
    MRE_number = 0
    pairs_key = []
    result_dict = {}
    mirna_length = len(mirna_seq) #miRNA_sequence_length
    targetrna_length = len(targetrna_revcomp_seq) #TargetRNA_sequence_length
    #print(mirna_seq)
    #if targetrna_length < 100:
    #    return result_dict
    mirna_seed = mirna_seq[0:seed_length] #miRNA_seed_sequence
    '''
    mirna_subseed_1 = ''
    mirna_subseed_2 = ''
    mirna_seed_rev = ''
    if len(mirna_seq) >= 21:
        mirna_subseed_1 = mirna_seq[12:16] #miRNA_sub-seed_sequence(13-16nt)
        mirna_subseed_2 = mirna_seq[16:21] #miRNA_sub-seed_sequence(15-21nt)
    elif len(mirna_seq) >= 16:
        mirna_subseed_1 = mirna_seq[12:16] #miRNA_sub-seed_sequence(13-16nt)
    if len(mirna_seq) >= 22:
        mirna_seed_rev = mirna_seq[14:22]
    else:
        mirna_seed_rev = mirna_seq[mirna_length-8:]
    '''

    for i in range(mirna_start_pairing-1,targetrna_length-mirna_length+1+1): #Start => End
        targetrna_seed = targetrna_revcomp_seq[i:i+seed_length]
        targetrna_seed_rev = ''
        targetrna_start = str(targetrna_length - (i + 0))
        targetrna_end = str(targetrna_length - (i + 7))
        #targetrna_sub1_start = str(targetrna_length - (i + 0))
        #targetrna_sub1_end = str(targetrna_length - (i + 4))
        #targetrna_sub2_start = str(targetrna_length - (i + 0))
        #targetrna_sub2_end = str(targetrna_length - (i + 5))
        #print (targetrna_start)
        key = mirna_id + '||' + targetrna_id + '||' + targetrna_end + '||' + targetrna_start
        #print (key)
        seed_match, seed_type, wobbles = find_seed_pairings(mirna_seed,targetrna_seed,seed_length,allowed_gu_wobbles,allowed_mismatches) #0 => No matched | 1=> matched
        #print (seed_match)
        if not seed_match == 0:
            targetrna_mirna_match = targetrna_revcomp_seq[i:i+mirna_length]
            '''
            subseed_1_desc = []
            subseed_2_desc = []
            targetrna_subseed_1 = ''
            targetrna_subseed_2 = ''
            #print(targetrna_end,targetrna_length)
            if(int(targetrna_end)+22) <= targetrna_length:
                #print('OK!')
                if mirna_subseed_1 != '':
                    targetrna_subseed_names1 = ['subseed1-2','subseed1-1','subseed1','subseed1+1','subseed1+2']
                    targetrna_subseed_1 = [targetrna_revcomp_seq[i+12-2:i+12+4-2],targetrna_revcomp_seq[i+12-1:i+12+4-1],targetrna_revcomp_seq[i+12:i+12+4],targetrna_revcomp_seq[i+12+1:i+12+4+1],targetrna_revcomp_seq[i+12+2:i+12+4+2]]
                    for x in range(0,len(targetrna_subseed_1)-1):
                        #print(targetrna_subseed_1[x],mirna_subseed_1)
                        if str(targetrna_subseed_1[x]) == mirna_subseed_1:
                            subseed_1_desc.append(targetrna_subseed_names1[x])
                    if len(subseed_1_desc) == 0:
                        subseed_1_desc.append('NA')
                else:
                    subseed_1_desc.append('NA')
                if mirna_subseed_2 != '':
                    targetrna_subseed_names2 = ['subseed2-2','subseed2-1','subseed2','subseed2+1','subseed2+2']
                    targetrna_subseed_2 = [targetrna_revcomp_seq[i+16-2:i+16+5-2],targetrna_revcomp_seq[i+16-1:i+16+5-1],targetrna_revcomp_seq[i+16:i+16+5],targetrna_revcomp_seq[i+16+1:i+16+5+1],targetrna_revcomp_seq[i+16+2:i+16+5+2]]
                    for x in range(0,len(targetrna_subseed_2)-1):
                        #print(targetrna_subseed_2[x],mirna_subseed_2)
                        if str(targetrna_subseed_2[x]) == mirna_subseed_2:
                            subseed_2_desc.append(targetrna_subseed_names2[x])
                    if len(subseed_2_desc) == 0:
                        subseed_2_desc.append('NA')
                else:
                    subseed_2_desc.append('NA')
            else:
                subseed_1_desc.append('NA')
                subseed_2_desc.append('NA')
            '''
            result_dict[key] = [mirna_seq, targetrna_mirna_match, mirna_seed,targetrna_seed,seed_match,seed_type,wobbles]
            pairs_key.append(key)
            MRE_number += 1
            #print (result_dict)
    if MRE_number == 0:
        return 0, 0
    site_density_simple_score = site_density_simple(MRE_number, targetrna_length)
    site_density_MuTaMe_score = site_density_MuTaMe(pairs_key, mirna_length)
    for key in pairs_key:
        result_dict[key].append(site_density_simple_score)
        result_dict[key].append(site_density_MuTaMe_score)

    return result_dict, pairs_key

if __name__ == '__main__':
    main()
