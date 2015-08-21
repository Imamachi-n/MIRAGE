#!usr/bin/env python

import os
import re
from parameter.common_parameters import common_parameters
import utils.setting_utils as utils

utils.now_time("phylop_score_prep script starting...")
p = utils.Bunch(common_parameters)

def main():
    utils.now_time("Input_file: " + p.phastcons_prep_input)
    utils.now_time("Output_file: " + p.phastcons_prep_output)

    for x in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']:
        input_s = p.phylop_score_prep_input + x + '.phyloP46way.wigFix'
        output_s = p.phylop_score_prep_input + x + '.phyloP46way.bed'
        phylop_score_prep_input_file = open(input_s,'r')
        phylop_score_prep_output_file = open(output_s,'w')

        chrom = ''
        start_site = 0
        step = 1

        for line in phylop_score_prep_input_file:
            line = line.rstrip()
            if re.match(r'^fixedStep',line):
                regex = r'fixedStep chrom=(?P<chrom>.+) start=(?P<start>.+) step=(?P<step>.+)'
                seq = re.match(regex,line)
                chrom = seq.group('chrom')
                start_site = int(seq.group('start')) - 1
                step = int(seq.group('step'))
                continue
            score = line
            #end_site = start_site + step
            for x in range(step):
                print (start_site, score, file=phylop_score_prep_output_file, sep="\t", end="\n")
                start_site += 1

        utils.now_time("phylop_score_prep script was successfully finished!!")
        phylop_score_prep_input_file.close()
        phylop_score_prep_output_file.close()

if __name__ == '__main__':
    main()
