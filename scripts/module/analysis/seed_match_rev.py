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

def is_gu_wobbles(mirna_nuc,targetrna_nuc):
    #check if mirna-targetrna pair is GU-wobble or not.
    if (mirna_nuc == 'G' and targetrna_nuc == 'A') or (mirna_nuc == 'U' and targetrna_nuc == 'C'):
        return True
    else:
        return False

def find_mirna_subtarget_candidates(mirna_id, mirna_seq, targetrna_id, targetrna_revcomp_seq, mirna_start_pairing=1, seed_length=8, allowed_gu_wobbles={6:0,7:1,8:2}, allowed_mismatches={6:0,7:0,8:0}):
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

    #Remove too short targetrna
    if targetrna_length < mirna_length:
        return 0, 0
    
    mirna_seed = mirna_seq[-8:] #miRNA_seed_sequence

    for i in range(mirna_start_pairing-1+mirna_length-8,targetrna_length-8): #Start => End
        targetrna_seed = targetrna_revcomp_seq[i:i+seed_length]
        targetrna_start = str(targetrna_length - (i + 0))
        targetrna_end = str(targetrna_length - (i + 7))
        key = mirna_id + '||' + targetrna_id + '||' + targetrna_end + '||' + targetrna_start
        seed_match = ''
        #print(mirna_seed,targetrna_seed)
        if mirna_seed == targetrna_seed:
            seed_match = '8mer_rev'
        else:
            seed_match = 0
        if not seed_match == 0:
            targetrna_mirna_match = targetrna_revcomp_seq[i:i+mirna_length]
            result_dict[key] = [mirna_seq, targetrna_mirna_match, mirna_seed,targetrna_seed,seed_match]
            pairs_key.append(key)
            MRE_number += 1
            #print (result_dict)
    if MRE_number == 0:
        return 0, 0
    return result_dict, pairs_key

if __name__ == '__main__':
    main()
