#!usr/bin/env python

import re
from parameter.common_parameters import common_parameters
import utils.setting_utils as utils

utils.now_time("Refseq_pre script starting...")
p = utils.Bunch(common_parameters)

def main():
    utils.now_time("Input_file: " + p.refseq_pre_input)
    utils.now_time("Output_file: " + p.refseq_pre_output)
    input_file = open(p.refseq_pre_input,'r')
    output_file = open(p.refseq_pre_output,'w')
    flg = 0
    seq = ""
    for line in input_file:
        line = line.rstrip()
        if re.match(r"^>",line): #Header
            data = line.split()
            refseq_id = data[0]
            refseq_id = refseq_id.replace('>hg19_refGene_','')
            if flg == 1:
                print (seq,file=output_file,end="\n")
            print (refseq_id,file=output_file,end="\t")
            flg = 1
            seq = ""
        else: #Sequence
            seq += line
    print (seq,file=output_file,end="\n")
    utils.now_time("Refseq_pre script was successfully finished!!")
    input_file.close()
    output_file.close()

if __name__ == '__main__':
    main()
