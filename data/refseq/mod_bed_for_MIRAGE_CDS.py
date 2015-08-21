#!/usr/bin/env python

import sys
import re

input_file = open('./refGene_2015-07-16_CDS.bed','r')
output_file = open('./refGene_2015-07-16_CDS_for_MIRAGE.bed','w')

checker = {}

for line in input_file:
    line = line.rstrip()
    data = line.split("\t")
    name = data[3]
    if not name in checker and re.match('^NM_',name):
        print(line, end="\n",file=output_file)
        checker[name] = 1
    else:
        continue

input_file.close()
output_file.close()
