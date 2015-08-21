#!/usr/bin/env python

from math import log2

ref_file = open('../../../result/mirage_output_miR-155_vs_RefSeq_NM_2015-07-17_miR-155_overexpression.result','r')
input_file = open('../../../result/gene_exp_miR-155_overexpression_RefSeq_Rep_isoforms.diff','r')
output_file = open('../../../result/gene_exp_miR-155_overexpression_RefSeq_Rep_isoforms_with_MRE_numbers_all_infor.diff','w')

ref_dict = {}
infor_dict = {}

for line in ref_file:
    line = line.rstrip()
    data = line.split("\t")
    gene_symbol = data[1]
    if data[0] == 'gr_id':
        infor_dict[gene_symbol] = "\t".join(data[14:])
        continue
    seed_type = data[22]
    GU_wobble = data[23]
    if GU_wobble == '0':
        if not gene_symbol in ref_dict:
            ref_dict[gene_symbol] = [seed_type]
            infor_dict[gene_symbol] = "\t".join(data[14:])
        else:
            ref_dict[gene_symbol].append(seed_type)
            infor_dict[gene_symbol] = 'NA'

for line in input_file:
    line = line.rstrip()
    data = line.split("\t")
    gene_symbol = data[1]
    if data[0] == 'gr_id':
        infor = infor_dict[gene_symbol]
        print(line,'seed_type','number',infor, sep="\t",end="\n",file=output_file)
        continue
    #print(data[0])
    value_1 = float(data[5])
    value_2 = float(data[6])
    if not gene_symbol in ref_dict:
        print(line,'NA',0, sep="\t",end="\n",file=output_file)
        continue
    else:
        seed_match_list = ref_dict[gene_symbol]
        infor = infor_dict[gene_symbol]
        number = len(seed_match_list)
        score = 0
        print(line,','.join(seed_match_list),number,infor, sep="\t",end="\n",file=output_file)

input_file.close()