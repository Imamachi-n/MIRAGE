#!/usr/bin/env python
'''
thermo_calc.py:
Calculate thermodynamic stability (minimum free energy (mfe) structures)


<Energy>
(1)miRNA seed region vs TargetRNA seed region
-------- miRNA(8nt_seed)
|||||||| 
-------- TargetRNA(8nt_seed)

(2)mature miRNA vs candidate target site (the same length)
---------------------- miRNA
||||||||||||||||||||||
---------------------- TargetRNA

(3)mature miRNA vs local TargetRNA region (70nt window)
        ---------------------- miRNA
        ||||||||||||||||||||||
-------------------------------------- TargetRNA

<Reference>
[1] Stormo GD. An overview of RNA structure prediction and applications to RNA gene prediction and RNAi design. Curr Protoc Bioinformatics. 2006 Mar;Chapter 12:Unit 12.1.
[2] http://www.tbi.univie.ac.at/RNA/tutorial/node6.html

'''

import shlex
import subprocess
import tempfile
import re

def make_constraints(seed_match, seq_type):
    if seq_type == 'miseq':
        seed_match = seed_match.replace('x','.')
        seed_match = seed_match.replace('A','.')
        seed_match = seed_match.replace(':','(')
        seed_match = seed_match.replace('|','(')
        return seed_match
    elif seq_type == 'targetseq': #Reverse
        seed_match_rev = seed_match[-1::-1]
        seed_match_rev = seed_match_rev.replace('x','.')
        seed_match_rev = seed_match_rev.replace('A','.')
        seed_match_rev = seed_match_rev.replace(':',')')
        seed_match_rev = seed_match_rev.replace('|',')')
        return  seed_match_rev

def viennaRNA_RNAcofold(seqs, constraints, option_postscript=False, option_constraints=True, option_partfunc=True, option_temperature=True):
    command_RNAcofold = 'RNAcofold --noPS --constraint --partfunc --temp=37'
    args = shlex.split(command_RNAcofold)
    test_str = "\n".join([seqs, constraints]) + "\n\n" + '@' + '\n'
    p = subprocess.Popen(args,stdin=subprocess.PIPE, stdout=subprocess.PIPE, cwd=tempfile.gettempdir())
    stdout, stderr = p.communicate("\n".join([seqs, constraints]).encode('utf-8')) #cannot use type 'str' for communicate...
    return stdout, stderr

#b'UUCAAGUA&UACUUGAA\n.(((((((&))))))). (-16.90)\n,(((((((&))))))), [-17.56]\n frequency of mfe structure in ensemble 0.342276 , delta G binding= -7.56\n'
def regex_RNAcofold(seq):
    regex = r'.+\n(?P<str_mfe>\S+) \((?P<mfe>.+)\)\n(?P<str_ens>\S+) \[(?P<ens>.+)\]\n frequency of mfe structure in ensemble (?P<ens_frequency>\S+) , delta G binding=(?P<delta_G>.+)\n'
    seq = seq.decode('utf-8')
    #print (seq)
    decoded_seq = re.match(regex, seq)
    str_mfe = decoded_seq.group('str_mfe')
    mfe = decoded_seq.group('mfe')
    str_ens = decoded_seq.group('str_ens')
    ens = decoded_seq.group('ens')
    delta_G = decoded_seq.group('delta_G')
    return str_mfe, mfe, str_ens, ens, delta_G

def calc_thermo(mirna_seq, targetrna_seq, targetrna_range, tmp_dict):
    mirna_length = len(mirna_seq) #miRNA sequence length
    targetrna_length = len(targetrna_seq)
    targetrna_range = 30 #Searched around 70nt
    around_nt_right = ''
    around_nt_left = ''
    if mirna_length % 2 == 0: #Even number
        around_nt_right = int((targetrna_range - mirna_length) / 2)
        around_nt_left = int((targetrna_range - mirna_length) / 2)
    else: #Odd number
        around_nt_right = int((targetrna_range - mirna_length - 1) / 2)
        around_nt_left = int((targetrna_range - mirna_length + 1) / 2)
    #miRNA_region
    mirna_seed = mirna_seq[0:8] #miRNA_seed_region
    mature_mirna = mirna_seq #mature_miRNA
    thermo_targetseq = '' #TargetRNA sequence for thermo calc.
    for x in list(tmp_dict.keys()):
        #print(x)
        mirna_infor = x
        mirna_data = mirna_infor.split('||')
        mirna_name = mirna_data[0]
        #TargetRNA_st_ed
        targetrna_ed = int(mirna_data[3]) #1nt - seed_region / end_site for miRNA-binding
        targetrna_st = targetrna_ed - mirna_length + 1 #8nt - seed_region / start_site for miRNA-binding
        #if (targetrna_st - around_nt_right) <= 0:
        #    print ('WARNINGS: ' + x)
        #    continue
        #if (targetrna_ed + around_nt_left) > targetrna_length:
        #    print ('WARNINGS: ' + x)
        #    continue
        #thermo_targetseq_st = targetrna_st - around_nt_right - 1
        #thermo_targetseq_ed = targetrna_ed + around_nt_left
        #Targetrna_region
        candidate_target_site = ''
        if not targetrna_st-1 < 0:
            candidate_target_site = targetrna_seq[targetrna_st-1:targetrna_ed]
        else:
            candidate_target_site = 'NA'
        #print(targetrna_st-1)
        targetrna_seed_region = targetrna_seq[targetrna_ed-8:targetrna_ed]
        #local_targetrna_region = targetrna_seq[thermo_targetseq_st:thermo_targetseq_ed] #TargetRNA sequence for thermo calc.
        #Calculated pairs
        test_seq1 = '&'.join([mirna_seed,targetrna_seed_region])
        test_seq2 = '&'.join([mature_mirna,candidate_target_site])
        #test_seq3 = '&'.join([mature_mirna,local_targetrna_region])

        #constraints
        c_miseq = ''
        c_targetseq = ''
        seed_match = (tmp_dict[x])[4] #NEED TO CHECK
        reside_miseq_targetseq = mirna_length - 8 #miseq - seed_region

        seed_match_miseq = make_constraints(seed_match,'miseq')
        c_miseq_seed = seed_match_miseq
        c_miseq = seed_match_miseq + reside_miseq_targetseq * '.'

        seed_match_targetseq = make_constraints(seed_match,'targetseq')
        c_targetseq_seed = seed_match_targetseq
        c_targetseq_site = reside_miseq_targetseq * '.' + seed_match_targetseq
        #c_targetseq = around_nt_right * '.' + reside_miseq_targetseq * '.' + seed_match_targetseq + around_nt_left * '.'

        test_constraints1 = '&'.join([c_miseq_seed,c_targetseq_seed])
        test_constraints2 = '&'.join([c_miseq,c_targetseq_site])
        #test_constraints3 = '&'.join([c_miseq,c_targetseq])

        #debug
        #print (test_seq1)
        #print (test_constraints1)
        #print (test_seq2)
        #print (test_constraints2)
        #print (test_seq3)
        #print (test_constraints3)

        #RNAcofold_command
        stdout1, stderr1 = viennaRNA_RNAcofold(test_seq1, test_constraints1) #test1
        stdout2 = ''
        stderr2 = ''
        if not candidate_target_site == 'NA':
            stdout2, stderr2 = viennaRNA_RNAcofold(test_seq2, test_constraints2) #test2
        else:
            stdout2 = 'NA'
        #stdout3, stderr3 = viennaRNA_RNAcofold(test_seq3, test_constraints3) #Test3
        #print (stdout1)
        #print (stdout2)
        #print (stdout3)
        #print (stderr)

        #Seed_matching
        str_mfe_seed, mfe_seed, str_ens_seed, ens_seed, delta_G_seed = regex_RNAcofold(stdout1)
        mfe_seed = mfe_seed.strip()
        ens_seed = ens_seed.strip()
        delta_G_seed = delta_G_seed.strip()
        out1_list = [str_mfe_seed, mfe_seed, str_ens_seed, ens_seed, delta_G_seed]
        tmp_dict[x].extend(out1_list)

        #miRNA-target_site matching
        if not stdout2 == 'NA':
            str_mfe, mfe, str_ens, ens, delta_G = regex_RNAcofold(stdout2)
            mfe = mfe.strip()
            ens = ens.strip()
            delta_G = delta_G.strip()
            out2_list = [str_mfe, mfe, str_ens, ens, delta_G]
            tmp_dict[x].extend(out2_list)

            #3'pairing contribution
            diff_mfe = float(mfe) - float(mfe_seed)
            diff_ens = float(ens) - float(ens_seed)
            diff_delta_G = float(delta_G) - float(delta_G_seed)
            out3_list = [diff_mfe, diff_ens, diff_delta_G]
            tmp_dict[x].extend(out3_list)
        else:
            tmp_dict[x].extend(['near_stop_codon','NA','NA','NA','NA'])
            tmp_dict[x].extend(['NA','NA','NA'])

        #print ('str_mfe: ' + str_mfe)
        #print ('mfe: ' + mfe)
        #print ('str_ens: ' + str_ens)
        #print ('ens: ' + ens)
        #print ('delta_G: ' + delta_G)
    return tmp_dict
