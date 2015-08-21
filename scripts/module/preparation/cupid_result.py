#!usr/bin/env python

import re
from parameter.common_parameters import common_parameters
from parameter.convert_mirbase_id import convert_mirbase_id
import utils.setting_utils as utils

utils.now_time("cupid_result script starting...")
p = utils.Bunch(common_parameters)

def main():
    utils.now_time("Input_file: " + p.cupid_pos)
    utils.now_time("Output_file: " + p.cupid_output)
    utils.now_time("miRNA_file: " + p.cupid_mirna_fasta)
    utils.now_time("targetRNA_file: " + p.cupid_targetrna_fasta)
    utils.now_time("Refseq_data: " + p.refseq_pre_output)
    utils.now_time("miRBase_data: " + p.mirbase_pre_output)
    refseq_dict = {}
    mirbase_dict = {}

    #mirbase_dict
    mirbase_file = open(p.mirbase_pre_output,'r')
    for line in mirbase_file:
        line = line.rstrip()
        data = line.split("\t")
        infor = data[0].split('|')
        mirbase_id = infor[0]
        symbol = infor[1]
        seq = data[1]
        if not re.match('hsa',mirbase_id):
            continue
        mirbase_dict[mirbase_id] = [symbol,seq] #miRNA_symbol => [0] | seq => [1]

    #refseq_dict 
    refseq_file = open(p.refseq_pre_output,'r')
    for line in refseq_file:
        line = line.rstrip()
        data = line.split("\t")
        refseq_id = data[0]
        seq = data[1]
        refseq_dict[refseq_id] = seq

    #main
    input_file = open(p.cupid_pos,'r')
    output_file = open(p.cupid_output,'w')
    mirna_file = open(p.cupid_mirna_fasta,'w')
    targetrna_file = open(p.cupid_targetrna_fasta,'w')
    error_file = open(p.cupid_error,'w')
    mirna_dist = {}
    targetrna_dist = {}
    for line in input_file:
        line = line.rstrip()
        data = line.split("\t")
        if data[0] == 'AvgProb[0,1]':
            continue
        mirbase_id = data[4]
        refseq_id = data[3]
        if refseq_id == "NM_000927":
            continue
        utr_infor = data[5].split('-')
        utr_st = int(utr_infor[0])
        utr_ed = int(utr_infor[1])
        if mirbase_id in convert_mirbase_id:
            mirbase_id = convert_mirbase_id[mirbase_id]
        if (refseq_id in refseq_dict and mirbase_id in mirbase_dict):
            symbol = mirbase_dict[mirbase_id][0]
            mir_seq = mirbase_dict[mirbase_id][1]
            mir_seq_length = len(mir_seq)
            utr_st = utr_ed - mir_seq_length - 5
            ref_seq_raw = refseq_dict[refseq_id]
            ref_seq = refseq_dict[refseq_id][utr_st:utr_ed]
            ref_seq = ref_seq.replace("T","U")
            mir_tag = '>' + mirbase_id + '|' + symbol
            refseq_tag = '>' + refseq_id
            print(mirbase_id,symbol,mir_seq,refseq_id,utr_st,utr_ed,ref_seq,file=output_file,sep="\t",end="\n")
            mirna_dist[mir_tag] = mir_seq
            targetrna_dist[refseq_tag] = ref_seq_raw
        else:
            print ("ERROR: " + refseq_id + '|' + mirbase_id,file=error_file,end="\n")
    for key in list(mirna_dist.keys()):
        print(key,file=mirna_file,end="\n")
        print(mirna_dist[key],file=mirna_file,end="\n")
    for key in list(targetrna_dist.keys()):
        print(key,file=targetrna_file,end="\n")
        print(targetrna_dist[key],file=targetrna_file,end="\n")

    utils.now_time("cupid_result script was successfully finished!!")
    input_file.close()
    output_file.close()

if __name__ == '__main__':
    main()
