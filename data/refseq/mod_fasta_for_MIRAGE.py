#!/usr/bin/env python

import sys
import re

input_file = './refGene_2015-06-03_3UTR.fa'
output_file = open('./refGene_2015-06-03_3UTR_for_MIRAGE.fasta','w')

#hg19_refGene_NM_032291 range=chr1:67208779-67210768 5'pad=0 3'pad=0 strand=+ repeatMasking=none
def load_fasta(fasta_path):
    fasta_file = open(fasta_path,'r')
    fasta_dict = {}
    checker = []
    name = ''
    flg = 0
    for line in fasta_file:
        line = line.rstrip()
        if line.startswith('#'):
            continue
        if line.startswith('>'):
            name = line[1:].strip()
            name_infor = name.split(" ")
            refid = name_infor[0].replace('hg19_refGene_','')
            if not refid in checker:
                checker.append(refid)
                fasta_dict[refid] = ''
                flg = 0
            else:
                #print ('ERROR: The same name exists =>' + refid)
                #sys.exit(1)
                flg = 1
        else:
            if flg == 1:
                continue
            trans_table = str.maketrans("ATGCatgcUu","AUGCAUGCUU")
            line_changed = line.translate(trans_table) #Translate DNA into RNA
            fasta_dict[refid] += line_changed.upper()
    fasta_file.close()
    return fasta_dict

fasta_dict = load_fasta(input_file)
for refid in list(fasta_dict.keys()):
    if re.match('^NR',refid):
        continue
    seq = fasta_dict[refid]
    print('>' + refid, sep="",end="\n",file=output_file)
    print(seq, end="\n",file=output_file)

output_file.close()
