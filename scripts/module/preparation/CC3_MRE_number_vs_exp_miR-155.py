#!/usr/bin/env python

from math import log2

ref_file = open('../../../result/mirage_output_miR-155_vs_RefSeq_NM_2015-07-17_miR-155_overexpression.result','r')
input_file = open('../../../result/gene_exp_miR-155_overexpression_RefSeq_Rep_isoforms.diff','r')
output_file = open('../../../result/gene_exp_miR-155_overexpression_RefSeq_Rep_isoforms_with_MRE_numbers.diff','w')
output_file2 = open('../../../result/gene_exp_miR-155_overexpression_RefSeq_Rep_isoforms_with_MRE_numbers_rpkm1.diff','w')

ref_dict = {}

mer6_m8 = 0.0929052-0.153277
mer6_m7 = -0.0578111-0.153277
mer7_m1 = -0.132557-0.153277
mer7_m8 = -0.236618-0.153277
mer8 = -0.37331-0.153277

for line in ref_file:
    line = line.rstrip()
    data = line.split("\t")
    if data[0] == 'gr_id':
        #print(line,end="\n",file=output_file)
        continue
    gene_symbol = data[1]
    seed_type = data[22]
    GU_wobble = data[23]
    if GU_wobble == '0':
        if not gene_symbol in ref_dict:
            ref_dict[gene_symbol] = [seed_type]
        else:
            ref_dict[gene_symbol].append(seed_type)

for line in input_file:
    line = line.rstrip()
    data = line.split("\t")
    if data[0] == 'gr_id':
        print(line,'seed_type','number','score', sep="\t",end="\n",file=output_file)
        print(line,'seed_type','number','score_log2', sep="\t",end="\n",file=output_file2)
        continue
    #print(data[0])
    gene_symbol = data[1]
    value_1 = float(data[5])
    value_2 = float(data[6])
    if not gene_symbol in ref_dict:
        print(line,'NA',0,0, sep="\t",end="\n",file=output_file)
        continue
    else:
        seed_match_list = ref_dict[gene_symbol]
        number = len(seed_match_list)
        score = 0
        for x in seed_match_list:
            if x == '8mer':
                score += mer8
            elif x == '7mer-m8':
                score += mer7_m8
            elif x == '7mer-m1':
                score += mer7_m1
            elif x == '6mer-m7':
                score += mer6_m7
            elif x == '6mer-m8':
                score += mer6_m8
        print(line,','.join(seed_match_list),number,score, sep="\t",end="\n",file=output_file)
        score_log2 = log2(-score)
        if value_1 >= 1 and value_2 >= 1:
            print(line,','.join(seed_match_list),number,score_log2, sep="\t",end="\n",file=output_file2)

input_file.close()