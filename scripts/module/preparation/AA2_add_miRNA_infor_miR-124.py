#!/usr/bin/env python

import re
from operator import itemgetter

ref_file = open('../../../data/RNA-seq_miR-124_miR-155_transfected_HeLa/gene_exp_miR-124_overexpression_RefSeq_Rep_isoforms.diff','r')
input_file = open('../../../result/mirage_output_miR-124_vs_RefSeq_NM_2015-07-20.result','r')
output_file = open('../../../result/mirage_output_miR-124_vs_RefSeq_NM_2015-07-20_miR-124_overexpression.result','w')

ref_dict = {}
header = ''

for line in ref_file:
    line = line.rstrip()
    data = line.split("\t")
    if data[0] == 'gr_id':
        header = line
        continue
    refid = data[2]
    ref_dict[refid] = line

for line in input_file:
    line = line.rstrip()
    data = line.split("\t")
    if data[0] == 'miRNA_name_id':
        print(header,line, sep="\t",end="\n",file=output_file)
        continue
    refid = data[1]
    if refid in ref_dict:
        print(ref_dict[refid],line, sep="\t",end="\n",file=output_file)

ref_file.close()
input_file.close()
output_file.close()
