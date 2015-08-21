#!/usr/bin/env python
'''
select_rep_isoforms_from_cuffdiff_output.py:
Select representative isoforms from cuffdiff output files

'''
import re
from operator import itemgetter

ref_file = open('../../../data/RNA-seq_miR-124_miR-155_transfected_HeLa/isoform_exp_miR-155_overexpression_RefSeq.diff','r')
input_file = open('../../../data/RNA-seq_miR-124_miR-155_transfected_HeLa/gene_exp_miR-155_overexpression_RefSeq.diff','r')
output_file = open('../../../data/RNA-seq_miR-124_miR-155_transfected_HeLa/gene_exp_miR-155_overexpression_RefSeq_Rep_isoforms.diff','w')

ref_dict = {}

for line in ref_file:
    line = line.rstrip()
    data = line.split("\t")
    if re.match('^test_id',data[0]):
        continue
    refid = data[0]
    gene_symbol = data[1]
    value_1 = data[7]
    if not gene_symbol in ref_dict:
        ref_dict[gene_symbol] = [[float(value_1), refid]]
    else:
        ref_dict[gene_symbol].append([float(value_1), refid])

gr_id = 1

for line in input_file:
    line = line.rstrip()
    data = line.split("\t")
    gene_symbol = data[0]
    locus = data[3]
    if data[0] == '':
        continue
    if re.match('^test_id',data[0]):
        print('gr_id','gene_symbol','refid','locus',"\t".join(data[6:]), sep="\t",end="\n",file=output_file)
        continue
    #print(ref_dict[gene_symbol])
    refid = sorted(ref_dict[gene_symbol],key=itemgetter(0),reverse=True)[0][1]
    if re.match('^NR_',refid):
        continue
    print(gr_id,gene_symbol,refid,locus,"\t".join(data[6:]), sep="\t",end="\n",file=output_file)
    gr_id += 1

ref_file.close()
input_file.close()
output_file.close()
