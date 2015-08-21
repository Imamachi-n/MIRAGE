#!usr/bin/env python

import sys
import re
import shelve
from parameter.common_parameters import common_parameters
import utils.setting_utils as utils

utils.now_time("phylop_score_list script starting...")
p = utils.Bunch(common_parameters)

def main():
    utils.now_time("Input_file: " + p.phylop_score_list_db_input)
    utils.now_time("Reference_file: " + p.phylop_score_list_reference)
    utils.now_time("Output_file: " + p.phylop_score_list_db_output)

    output_merge = p.phylop_score_list_db_output + 'phyloP46way_Refseq_for_MIRAGE_CDS.db' #'phyloP46way_miRBase_v21_hg38Tohg19_for_MIRAGE.db'
    output_merge_shelve = shelve.open(output_merge)

    #for x in ['chrY']:
    for x in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']:
        ref_s = p.phylop_score_list_reference #mirBase, Refseq etc...
        ref_file = open(ref_s,'r')

        input_s = p.phylop_score_list_db_input + x + '.phyloP46way_Refseq_CDS.db' #'.phyloP46way_Refseq.db'
        output_s = p.phylop_score_list_db_output +  x + '.phyloP46way_Refseq_for_MIRAGE_CDS.db' #'.phyloP46way_Refseq_for_MIRAGE.db'

        input_shelve = shelve.open(input_s)
        output_shelve = shelve.open(output_s)

        score_list_dict = {}

        for line in ref_file:
            line = line.rstrip()
            data = line.split("\t")
            chrom = data[0]
            if not chrom == x:
                continue
            strand = data[5]
            if len(data) >= 12: #12bed format
                exon_block = data[10].split(',')
                exon_block.pop() #Remove the last item ''
                exon_st = data[11].split(',')
                exon_st.pop() #Remove the last item ''
                name = data[3]
                score_list_dict[name] = []
                for y in range(len(exon_block)):
                    st = int(data[1]) + int(exon_st[y])
                    ed = int(data[1]) + int(exon_st[y]) + int(exon_block[y])
                    length = ed - st
                    for z in range(length):
                        score = input_shelve[str(st)]
                        score_list_dict[name].append(score)
                        st += 1
                if strand == '-':
                    rev_score = score_list_dict[name][::-1]
                    score_list_dict[name] = rev_score
            elif len(data) >= 3: #6bed format
                st = int(data[1])
                ed = int(data[2])
                length = ed - st
                name = data[3]
                score_list_dict[name] = []
                for z in range(length):
                    score = input_shelve[str(st)]
                    score_list_dict[name].append(score)
                    st += 1
                if strand == '-':
                    rev_score = score_list_dict[name][::-1]
                    score_list_dict[name] = rev_score
            else:
                print('ERROR: Your BED format file have less than three column.')
                print ('BED format file need to have at least three column [chr, st, ed]...')
                sys.exit(1)

        output_shelve.update(score_list_dict)
        output_merge_shelve.update(score_list_dict)
        input_shelve.close()
        output_shelve.close()

    utils.now_time("phylop_score_list script was successfully finished!!")
    output_merge_shelve.close()

if __name__ == '__main__':
    main()
