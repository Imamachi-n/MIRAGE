#Title: phastcons_score_heatmap
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-06-10

#library
library(ggplot2)

#file_path
setwd("C:/Users/Naoto/Documents/github/MIRAGE/result/")
input_file_path <- "C:/Users/Naoto/Documents/github/MIRAGE/result/gene_exp_miR-124_overexpression_RefSeq_Rep_isoforms_with_MRE_numbers_rpkm1.diff"

input_file <- read.table(input_file_path,header=T)
colnames(input_file) <- c('gr_id','gene_symbol','refid','locus','status','value_1','value_2','log2_fold_change','test_stat',
'p_value','q_value','significant','seed_type','number','score_log2')

p.scatter <- ggplot(input_file,aes(x=log2_fold_change, y=score_log2))
p.scatter <- p.scatter +  geom_point(shape = 20,               # �v���b�g�̃^�C�v���w��
                                     size = 4.0,               # �v���b�g�̃T�C�Y���w��
                                     na.rm = TRUE              # �񐔒l�𖳎�
                                     )
p.scatter <- p.scatter + geom_smooth(       # �ߎ���
                                     method = "lm"             # �ߎ����͉�A�@�ɂ���ċ��߂�
                                    )
p.scatter <- p.scatter + xlab("log2(Fold_change)")    # x �����x��
p.scatter <- p.scatter + ylab("log2(MRE_score)")    # y �����x��
p.scatter <- p.scatter + ggtitle("MRE_score") # �O���t�^�C�g��

plot(p.scatter)