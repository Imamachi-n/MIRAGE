#!/usr/bin/env python

common_parameters = dict(
#bed_3UTR
bed_3UTR_input = '../data/refseq/refGene_2015-06-03.bed',
bed_3UTR_output = '../data/refseq/refGene_2015-06-03_3UTR.bed',

#MIRBASE_GFF2BED
mirbase_gff2bed_input = '../data/miRBase/v21/hsa.gff3',
mirbase_gff2bed_output = '../data/miRBase/v21/hsa.bed',

#liftOver
liftover_input = '../data/miRBase/v21/hsa.bed',
liftover_output = '../data/miRBase/v21/hsa_hg38Tohg19.bed',

#phylop_score_prep
phylop_score_prep_input = '../data/PhyloP/',
phylop_score_prep_output = '../data/PhyloP/',

#phastcons_prep
phastcons_prep_input = '../data/PhastCons46Ways/', #chr4.phastCons46way.wigFix
phastcons_prep_output = '../data/PhastCons46Ways/', #chr4.phastCons46way.bed

#phastcons_sizedown
#phastcons_sizedown_bed_input = '../data/refseq/refGene_2015-07-16_3UTR_for_MIRAGE.bed',
phastcons_sizedown_bed_input = '../data/refseq/refGene_2015-07-16_CDS.bed',
#phastcons_sizedown_bed_input = '../data/miRBase/v21/hsa_hg38Tohg19.bed',
phastcons_sizedown_score_input = '../data/PhastCons46Ways/',
phastcons_sizedown_score_output = '../data/PhastCons46Ways/',

#phylop_sizedown
#phylop_sizedown_bed_input = '../data/refseq/refGene_2015-07-16_3UTR_for_MIRAGE.bed',
phylop_sizedown_bed_input = '../data/refseq/refGene_2015-07-16_CDS.bed',
#phylop_sizedown_bed_input = '../data/miRBase/v21/hsa_hg38Tohg19.bed',
phylop_sizedown_score_input = '../data/PhyloP/',
phylop_sizedown_score_output = '../data/PhyloP/',

#phylop_score_list
#phylop_score_list_reference = '../data/refseq/refGene_2015-07-16_3UTR_for_MIRAGE.bed',
phylop_score_list_reference = '../data/refseq/refGene_2015-07-16_CDS.bed',
#phylop_score_list_reference = '../data/miRBase/v21/hsa_hg38Tohg19.bed',
phylop_score_list_db_input = '../data/PhyloP/',
phylop_score_list_db_output = '../data/PhyloP/',

#phastcons_score_list
#phastcons_score_list_reference = '../data/refseq/refGene_2015-07-16_3UTR_for_MIRAGE.bed',
phastcons_score_list_reference = '../data/refseq/refGene_2015-07-16_CDS.bed',
#phastcons_score_list_reference = '../data/miRBase/v21/hsa_hg38Tohg19.bed',
phastcons_score_list_db_input = '../data/PhastCons46Ways/',
phastcons_score_list_db_output = '../data/PhastCons46Ways/',

#phastcons_score_R
phastcons_score_R_input = '../data/PhastCons46Ways/',
phastcons_score_R_output = '../data/PhastCons46Ways/',

#phylop_score_R
phylop_score_R_input = '../data/PhyloP/',
phylop_score_R_output = '../data/PhyloP/',

#REFSEQ_PRE
refseq_pre_input = '../data/refseq/refGene_2014-04-23_3UTR.fasta',
refseq_pre_output = '../data/refseq/refGene_2014-04-23_3UTR_id-seq.txt',

#MIRBASE_PRE
mirbase_pre_input = '../data/miRBase/v21/mature.fa',
mirbase_pre_output = '../data/miRBase/v21/miRBase_v21_symbol-id-seq.txt',

#MIRMARK_POS
mirmark_pos = '../data/mirMark/S1_site-level_human_miRNA-target_site_pairs_from_miRecords.csv',
mirmark_output = '../data/mirMark/S1_site-level_human_miRNA-target_site_pairs_from_miRecords_id-seq.txt',
mirmark_error = '../data/mirMark/S1_site-level_human_miRNA-target_site_pairs_error.txt',
mirmark_mirna_fasta = '../data/mirMark/mirmark_mirna.fa',
mirmark_targetrna_fasta = '../data/mirMark/mirmark_targetrna.fa',

#CUPID_POS
cupid_pos = '../data/cupid/tablePos_1481sites.txt',
cupid_output = '../data/cupid/tablePos_1481sites_id-seq.txt',
cupid_error = '../data/cupid/tablePos_1481sites_error.txt',
cupid_mirna_fasta = '../data/cupid/cupid_mirna.fa',
cupid_targetrna_fasta = '../data/cupid/cupid_targetrna.fa',
)
