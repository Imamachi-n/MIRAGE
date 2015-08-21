#!/usr/bin/env python3

input_file = open('../../../data/PAR-CLIP_HuR/SRR309285_4SU_DMEM_PAR-CLIP_HuR_clusters.csv','r')
output_file = open('../../../data/PAR-CLIP_HuR/SRR309285_4SU_DMEM_PAR-CLIP_HuR_clusters.bed','w')

for line in input_file:
    line = line.rstrip()
    data = line.split(',')
    if data[0] == 'Chromosome':
        continue
    chrom = data[0]
    strand = data[1]
    st = data[2]
    ed = data[3]
    name = data[4]
    score = data[6]
    print(chrom,st,ed,name,score,strand, sep="\t",end="\n",file=output_file)

input_file.close()
output_file.close()


input_file = open('../../../data/PAR-CLIP_HuR/SRR309286_4SU_SILAC_PAR-CLIP_HuR_clusters.csv','r')
output_file = open('../../../data/PAR-CLIP_HuR/SRR309286_4SU_SILAC_PAR-CLIP_HuR_clusters.bed','w')

for line in input_file:
    line = line.rstrip()
    data = line.split(',')
    if data[0] == 'Chromosome':
        continue
    chrom = data[0]
    strand = data[1]
    st = data[2]
    ed = data[3]
    name = data[4]
    score = data[6]
    print(chrom,st,ed,name,score,strand, sep="\t",end="\n",file=output_file)

input_file.close()
output_file.close()
