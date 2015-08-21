#!/usr/bin/env python

ref_file = open('../../../result/PhastCons46Ways_miRBase_v21_hg38Tohg19_kmeans5_2015-07-17.txt','r')
input_file = open('../../../data/mirMark/S1_site-level_human_miRNA-target_site_pairs_from_miRecords_id-seq.txt','r')
output_file = open('../../../result/PhastCons46Ways_miRBase_v21_hg38Tohg19_with_miRecords_kmeans5_2015-07-20.txt','w')

ref_dict = {}
ref_counter = {}
total = 0

for line in ref_file:
    line = line.rstrip()
    data = line.split("\t")
    if data[2] == '2':
        continue
    mirna_infor = data[0].split('|')
    mirna_id = mirna_infor[1]
    score = data[1]
    if not score in ref_counter:
        ref_counter[score] = 1
        total += 1
    else:
        ref_counter[score] += 1
        total += 1
    if not mirna_id in ref_dict:
        ref_dict[mirna_id] = score
    else:
        continue

print('1: ',ref_counter['1'],int(ref_counter['1'])/total*100)
print('2: ',ref_counter['2'],int(ref_counter['2'])/total*100)
print('3: ',ref_counter['3'],int(ref_counter['3'])/total*100)
print('4: ',ref_counter['4'],int(ref_counter['4'])/total*100)
print('5: ',ref_counter['5'],int(ref_counter['5'])/total*100)

counter = {}
kmeans_count = {}
kmeans_count['3'] = 0
kmeans_count['4'] = 0

for line in input_file:
    line = line.rstrip()
    data = line.split("\t")
    mirna_id = data[1]
    score = ref_dict[mirna_id]
    if not mirna_id in counter:
        if not score in kmeans_count:
            kmeans_count[score] = 1
        else:
            kmeans_count[score] += 1
        print(mirna_id,score, sep="\t",end="\n",file=output_file)
        counter[mirna_id] = 1

print('1: ',kmeans_count['1'],int(kmeans_count['1'])/129*100)
print('2: ',kmeans_count['2'],int(kmeans_count['2'])/129*100)
print('3: ',kmeans_count['3'],int(kmeans_count['3'])/129*100)
print('4: ',kmeans_count['4'],int(kmeans_count['4'])/129*100)
print('5: ',kmeans_count['5'],int(kmeans_count['5'])/129*100)
