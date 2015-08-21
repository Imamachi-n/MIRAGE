#!/usr/bin/env python

import shelve
from parameter.common_parameters import common_parameters
import utils.setting_utils as utils

utils.now_time("phastcons_score_R script starting...")
p = utils.Bunch(common_parameters)

def main():
    utils.now_time("Input_file: " + p.phastcons_score_R_input)
    utils.now_time("Output_file: " + p.phastcons_score_R_output)

    output_s = p.phastcons_score_R_output + 'phastCons46way_miRBase_v21_hg38Tohg19.txt'
    output_file = open(output_s,'w')

    #for x in ['chrY']:
    for x in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']:
        input_s = p.phastcons_score_R_input + x + '.phastCons46way_miRBase_v21_hg38Tohg19_for_MIRAGE.db'
        input_shelve = shelve.open(input_s)
        max_length = 28 #Max_length: 28nt(miRNA)
        for keys in input_shelve.keys():
            values = input_shelve[keys]
            value_length = len(values)
            add_length = max_length - value_length
            null_value = [0.000 for i in range(add_length)]
            values += null_value
            value_string = "\t".join(map(str, values))
            print(keys,value_string, file=output_file, sep="\t", end="\n")
        input_shelve.close()

    output_file.close()

    utils.now_time("phastcons_score_R script was successfully finished!!")

if __name__ == '__main__':
    main()
