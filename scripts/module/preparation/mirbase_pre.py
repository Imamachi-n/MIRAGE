#!usr/bin/env python

import re
from parameter.common_parameters import common_parameters
import utils.setting_utils as utils

utils.now_time("mirbase_pre script starting...")
p = utils.Bunch(common_parameters)

def main():
    utils.now_time("Input_file: " + p.mirbase_pre_input)
    utils.now_time("Output_file: " + p.mirbase_pre_output)
    input_file = open(p.mirbase_pre_input,'r')
    output_file = open(p.mirbase_pre_output,'w')
    flg = 0
    seq = ""
    for line in input_file:
        line = line.rstrip()
        if re.match(r"^>",line): #Header
            data = line.split()
            mir_id = data[0]
            mir_id = mir_id.replace('>','')
            symbol = data[1]
            infor = mir_id + '|' + symbol
            if flg == 1:
                print (seq,file=output_file,end="\n")
            print (infor,file=output_file,end="\t")
            flg = 1
            seq = ""
        else: #Sequence
            seq += line
    print (seq,file=output_file,end="\n")
    utils.now_time("mirbase_pre script was successfully finished!!")
    input_file.close()
    output_file.close()

if __name__ == '__main__':
    main()
