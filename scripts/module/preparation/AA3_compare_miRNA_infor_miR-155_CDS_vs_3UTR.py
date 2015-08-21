#!/usr/bin/env python

import re
from operator import itemgetter

ref_file = open('../../../result/mirage_output_miR-155_vs_RefSeq_NM_2015-07-17_miR-155_overexpression.result','r')
ref_file2 = open('../../../result/genes_mRNA-seq.fpkm_table_merged_Ribo-seq.result','r')
input_file = open('../../../result/mirage_output_miR-155_vs_RefSeq_NM_CDS_2015-08-06_miR-155_overexpression.result','r')
output_file = open('../../../result/mirage_output_miR-155_vs_RefSeq_NM_CDS_2015-08-06_miR-155_overexpression_non_MRE_in_3UTR.result','w')

ref_dict = {}

for line in ref_file:
    line =line.rstrip()
    data = line.split("\t")
    if data[0] == 'gr_id':
        continue
    gr_id = data[0]
    ref_dict[gr_id] = 1

ref_dict2 = {}

for line in ref_file2:
    line = line.rstrip()
    data = line.split("\t")
    symbol = data[1]
    ref_dict2[symbol] = line

for line in input_file:
    line = line.rstrip()
    data = line.split("\t")
    if data[0] == 'gr_id':
        print("\t".join(data[0:26]),ref_dict2['Gene_ID'], sep="\t",end="\n",file=output_file)
        continue
    gr_id = data[0]
    symbol = data[1]
    if not gr_id in ref_dict:
        if symbol in ref_dict2:
            print("\t".join(data[0:26]),ref_dict2[symbol], sep="\t",end="\n",file=output_file)

ref_file.close()
ref_file2.close()
input_file.close()
output_file.close()
