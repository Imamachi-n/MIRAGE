#!usr/bin/env python

import sys
import re
import shelve
from parameter.common_parameters import common_parameters
import utils.setting_utils as utils

utils.now_time("phylop_sizedown script starting...")
p = utils.Bunch(common_parameters)

def main():
    utils.now_time("Input_file: " + p.phylop_sizedown_score_input)
    utils.now_time("Reference_file: " + p.phylop_sizedown_bed_input)
    utils.now_time("Output_file: " + p.phylop_sizedown_score_output)

    '''
    ref_s = p.phastcons_sizedown_bed_input #mirBase, Refseq etc...
    ref_file = open(ref_s,'r')
    ref_dict = {} #{NM_000XXXX: [st1,ed1],[st2,ed2]}
    for line in ref_file:
        line = line.rstrip()
        data = line.split("\t")

        if len(data) >= 12: #12bed format
            st = 0
            ed = 0
            exon_block = data[10].split(',')
            exon_block.pop()
            exon_st = data[11].split(',')
            exon_st.pop()
            chrom = data[0]
            name = data[3]
            for y in range(len(exon_block)):
                st = int(data[1]) + int(exon_st[y])
                ed = int(data[1]) + int(exon_st[y]) + int(exon_block[y])
                if not name in ref_dict:
                    ref_dict[name] = [[chrom,st,ed]]
                else:
                    ref_dict[name].append([chrom,st,ed])
        else: #6bed format
            st = data[1]
            ed = data[2]
            name = data[3]
            if not name in ref_dict:
                ref_dict[name] = [[chrom,st,ed]]
            else:
                ref_dict[name].append([chrom,st,ed])
    '''

    for x in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']:
    #for x in ['chrY']:
        ref_s = p.phylop_sizedown_bed_input #mirBase, Refseq etc...
        ref_file = open(ref_s,'r')

        input_s = p.phylop_sizedown_score_input + x + '.phyloP46way.bed'
        output_s = p.phylop_sizedown_score_output +  x + '.phyloP46way_Refseq_CDS.db'
        phylop_sizedown_input_file = open(input_s,'r')

        score_dict = {}

        for line in ref_file:
            line = line.rstrip()
            data = line.split("\t")
            chrom = data[0]
            if not x == chrom:
                continue
            if len(data) >= 12: #12bed format
                exon_block = data[10].split(',')
                exon_block.pop() #Remove the last item ''
                exon_st = data[11].split(',')
                exon_st.pop() #Remove the last item ''
                #name = data[3]
                for y in range(len(exon_block)):
                    st = int(data[1]) + int(exon_st[y])
                    ed = int(data[1]) + int(exon_st[y]) + int(exon_block[y])
                    length = ed - st
                    for z in range(length):
                        score_dict[str(st)] = 0
                        st += 1
            elif len(data) >= 3: #6bed format
                st = int(data[1])
                ed = int(data[2])
                length = ed - st
                for z in range(length):
                    score_dict[str(st)] = 0
                    st += 1
            else:
                print('ERROR: Your BED format file have less than three column.')
                print ('BED format file need to have at least three column [chr, st, ed]...')
                sys.exit(1)

        utils.now_time('Reference_file was loaded.')

        for line in phylop_sizedown_input_file:
            line = line.rstrip()
            data = line.split("\t")
            st_site = 0
            score = 0
            if re.match(r'^chr',data[0]):
                st_site = data[1] #
                score = data[2] #
            else:
                st_site = data[0] #
                score = data[1] #
            if st_site in score_dict:
                score_dict[str(st_site)] = score

        shelve_db = shelve.open(output_s)
        shelve_db.update(score_dict)
            
        utils.now_time("phylop_sizedown script was successfully finished!!")
        phylop_sizedown_input_file.close()
        shelve_db.close()


if __name__ == '__main__':
    main()
