#!/usr/bin/env python

"""
Compare_mirMark_and_seed_matching_result.

Created by Naoto Imamachi on 2015-03-22.

Usage:
  compare_mirMark_and_seed_matching_result.py <analysis_type> <miRNA.fasta> <targetRNA.fasta> [option]

"""

import argparse
import shelve

def is_gu_wobbles(mirna_nuc,targetrna_nuc):
    #check if mirna-targetrna pair is GU-wobble or not.
    if (mirna_nuc == 'G' and targetrna_nuc == 'U') or (mirna_nuc == 'U' and targetrna_nuc == 'G'):
        return True
    else:
        return False

def reverse_complement(seq):
    if 'U' in seq or 'u' in seq:
        trans_table = str.maketrans("AUGCaugc","UACGUACG")
    else:
        trans_table = str.maketrans("ATGCatgc","UACGUACG")
    rev_seq = seq[-1::-1]
    rev_comp_seq = rev_seq.translate(trans_table)
    return rev_comp_seq

def main():
    parser = argparse.ArgumentParser(prog='Compare_mirMark_and_seed_matching_result',description='Compare_experimental_verified_data_and_predicted_seed_matching_result')
    parser.add_argument('-i', '--input-result-file', action='store', required=True, dest='input_file', help='Specify your input result file')
    parser.add_argument('-r','-reference-verified-file', action='store', required=True, dest='verified_file', help='Specify your reference verified data [seed_matching_result.db]')
    parser.add_argument('-o','--output-file', action='store', required=True, dest='result_file', help='Sepecify your result data')

    args = parser.parse_args()
    input_file_path = args.input_file
    ref_file_path = args.verified_file
    output_file_path = args.result_file

    #'hsa-miR-1285-3p|MIMAT0005876||369||375': ['UCUGGGCA', 'CCUGAGCA', ':|||:|||', '8mer', 2]
    ref_file_data = open(ref_file_path,'r')
    ref_file = {}
    for x in ref_file_data:
        x.rstrip()
        data = x.split("\t")
        key = data[0] + '||' + data[1] + '||' + data[2] + '||' + data[3]
        value = data[4:]
        ref_file[key] = value
    input_file = open(input_file_path,'r')
    output_file = open(output_file_path,'w')
    for line in input_file:
        line = line.rstrip()
        data = line.split("\t")
        mir_symbol = data[0]
        mir_id = data[1]
        ref_id = data[3]
        #print (data[5])
        ref_st = int(data[5]) - 7
        ref_ed = int(data[5]) - 1
        key_name = mir_symbol + '|' + mir_id + '||' + ref_id + '||' + str(ref_st) + '||' + str(ref_ed)
        print (line,file=output_file,end="\t")
        if int(data[4]) < 0:
            print ("\t",file=output_file,end="")
        if ref_file.get(key_name,'None') == 'None': #Not found
            if int(data[4]) < 0:
                print ('NA' + "\t" + 'NA',file=output_file,end="\n")
                continue
            else:
                seq_mirna = data[2]
                seq_targetrna = reverse_complement(data[6])
                seq_mirna_seed = seq_mirna[0:8]
                seq_targetrna_seed = seq_targetrna[0:8]
                seq_type = ''
                for x in range(8):
                    if seq_mirna_seed[x] == seq_targetrna_seed[x]:
                        seq_type += '|'
                    elif is_gu_wobbles(seq_mirna_seed[x],seq_targetrna_seed[x]):
                        seq_type += ':'
                    else:
                        seq_type += 'x'
                print ('NA' + "\t" + 'NA' + "\t" + seq_type,file=output_file,end="\n")
        else: #Found
            value = ref_file[key_name]
            v_type = value[2]
            v_mer = value[3]
            if ':' in v_type: #Found by MIRAGE
                print ('NA' + "\t" + 'MIRAGE' + "\t" + v_type + "\t" + v_mer,file=output_file,end="\n")
            else: #Found by targetscan
                print ('TargetScan' + "\t" + 'MIRAGE' + "\t" + v_type + "\t" + v_mer,file=output_file,end="\n")

if __name__ == '__main__':
    main()
