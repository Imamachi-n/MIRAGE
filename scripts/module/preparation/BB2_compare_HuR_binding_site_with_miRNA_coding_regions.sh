#!/bin/bash
#HeLaS3
cd ../../..
bedtools intersect -b ./data/miRBase/v21/hsa_hg38Tohg19.bed -a ./data/PAR-CLIP_HuR/SRR309285_4SU_DMEM_PAR-CLIP_HuR_clusters.bed -wa -wb > ./data/PAR-CLIP_HuR/SRR309285_4SU_DMEM_PAR-CLIP_HuR_clusters_with_miRBase_v21.bed
bedtools intersect -b ./data/miRBase/v21/hsa_hg38Tohg19.bed -a ./data/PAR-CLIP_HuR/SRR309286_4SU_SILAC_PAR-CLIP_HuR_clusters.bed -wa -wb > ./data/PAR-CLIP_HuR/SRR309286_4SU_SILAC_PAR-CLIP_HuR_clusters_with_miRBase_v21.bed
