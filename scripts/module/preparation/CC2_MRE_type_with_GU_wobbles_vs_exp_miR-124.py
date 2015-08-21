#!/usr/bin/env python

ref_file = open('../../../result/mirage_output_miR-124_vs_RefSeq_NM_2015-07-20_miR-124_overexpression.result','r')
input_file = open('../../../result/gene_exp_miR-124_overexpression_RefSeq_Rep_isoforms.diff','r')
output_file = open('../../../result/gene_exp_miR-124_overexpression_RefSeq_Rep_isoforms_with_seed_type_wobbles.diff','w')

ref_dict = {}

for line in ref_file:
    line = line.rstrip()
    data = line.split("\t")
    if data[0] == 'gr_id':
        #print(line,end="\n",file=output_file)
        continue
    gene_symbol = data[1]
    seed_type = data[22]
    GU_wobble = data[23]
    seed_type2 = seed_type + '|' + GU_wobble
    if not gene_symbol in ref_dict:
        ref_dict[gene_symbol] = [seed_type2]
    else:
        ref_dict[gene_symbol].append(seed_type2)

#print(ref_dict['ZSWIM2'])

for line in input_file:
    line = line.rstrip()
    data = line.split("\t")
    if data[0] == 'gr_id':
        print(line,'seed_type', sep="\t",end="\n",file=output_file)
        continue
    #print(data[0])
    gene_symbol = data[1]
    if not gene_symbol in ref_dict:
        print(line,'NA', sep="\t",end="\n",file=output_file)
        continue
    else:
        seed_match_list = ref_dict[gene_symbol]
        if '8mer|0' in seed_match_list:
            print(line,'8mer|0', sep="\t",end="\n",file=output_file)
        elif '8mer|1' in seed_match_list:
            print(line,'8mer|1', sep="\t",end="\n",file=output_file)
        elif '8mer|2' in seed_match_list:
            print(line,'8mer|2', sep="\t",end="\n",file=output_file)
        elif '7mer-m8|0' in seed_match_list:
            print(line,'7mer-m8|0', sep="\t",end="\n",file=output_file)
        elif '7mer-m1|0' in seed_match_list:
            print(line,'7mer-m1|0', sep="\t",end="\n",file=output_file)
        elif '7mer-m8|1' in seed_match_list:
            print(line,'7mer-m8|1', sep="\t",end="\n",file=output_file)
        elif '7mer-m1|1' in seed_match_list:
            print(line,'7mer-m1|1', sep="\t",end="\n",file=output_file)
        elif '6mer-m7|0' in seed_match_list:
            print(line,'6mer-m7|0', sep="\t",end="\n",file=output_file)
        elif '6mer-m8|0' in seed_match_list:
            print(line,'6mer-m8|0', sep="\t",end="\n",file=output_file)
        else:
            #pass
            print('ERROR: ' + gene_symbol)

input_file.close()
