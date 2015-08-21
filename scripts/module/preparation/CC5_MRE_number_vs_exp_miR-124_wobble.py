#!/usr/bin/env python

from math import log2

ref_file = open('../../../result/mirage_output_miR-124_vs_RefSeq_NM_2015-07-20_miR-124_overexpression.result','r')
input_file = open('../../../result/gene_exp_miR-124_overexpression_RefSeq_Rep_isoforms.diff','r')
output_file = open('../../../result/gene_exp_miR-124_overexpression_RefSeq_Rep_isoforms_with_MRE_numbers_wobbles.diff','w')

ref_dict = {}


for line in ref_file:
    line = line.rstrip()
    data = line.split("\t")
    if data[0] == 'gr_id':
        #print(line,end="\n",file=output_file)
        continue
    gene_symbol = data[1]
    seed_type = data[22]
    GU_wobble = data[23]
    seed_type2 = seed_type + '|' + GU_wobble
    if not gene_symbol in ref_dict:
        ref_dict[gene_symbol] = [seed_type2]
    else:
        ref_dict[gene_symbol].append(seed_type2)

for line in input_file:
    line = line.rstrip()
    data = line.split("\t")
    if data[0] == 'gr_id':
        print(line,'seed_type','number', sep="\t",end="\n",file=output_file)
        continue
    #print(data[0])
    gene_symbol = data[1]
    if not gene_symbol in ref_dict:
        print(line,'NA',0, sep="\t",end="\n",file=output_file)
        continue
    else:
        seed_match_list = ref_dict[gene_symbol]
        number = len(seed_match_list)
        print(line,','.join(seed_match_list),number, sep="\t",end="\n",file=output_file)

input_file.close()
