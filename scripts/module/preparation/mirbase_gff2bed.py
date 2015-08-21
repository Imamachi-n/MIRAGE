#!usr/bin/env python

import re
from parameter.common_parameters import common_parameters
import utils.setting_utils as utils

utils.now_time("mirbase_gff2bed script starting...")
p = utils.Bunch(common_parameters)

def main():
    utils.now_time("Input_file: " + p.mirbase_gff2bed_input)
    utils.now_time("Output_file: " + p.mirbase_gff2bed_output)

    mirbase_gff_file = open(p.mirbase_gff2bed_input,'r')
    mirbase_bed_file = open(p.mirbase_gff2bed_output,'w')

    for line in mirbase_gff_file:
        line = line.rstrip()
        data = line.split("\t")
        if re.match(r'^#',line):
            continue
        chrom = data[0]
        status = data[2]
        st = int(data[3]) - 1
        ed = data[4]
        strand = data[6]
        if status == 'miRNA_primary_transcript':
            continue
        name_infor = data[8].split(';')
        mir_id = re.sub(r'^ID=','',name_infor[0])
        mir_id_number = ''
        if re.search(r'_',mir_id):
            mir_id, mir_id_number = mir_id.split('_')
        else:
            mir_id_number = 0 #there is ONLY one miRNA coding site in your genome
        mir_name = re.sub(r'^Name=','',name_infor[2])
        name = mir_name + '|' + mir_id + '|' + str(mir_id_number)
        print (chrom, st, ed, name, 0, strand, file=mirbase_bed_file, sep="\t", end="\n")

    utils.now_time("mirbase_gff2bed script was successfully finished!!")
    mirbase_gff_file.close()
    mirbase_bed_file.close()

if __name__ == '__main__':
    main()
