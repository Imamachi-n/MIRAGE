#phylop_score_variation_1_22_nt.r:
#Calculate average conservation score(PhyloP) in each nucleotide of miRNA(1-22nt)

#library
library(reshape2)
library(ggplot2)

data <- read.table('C:/Users/Naoto/Documents/github/MIRAGE/data/PhyloP/phyloP46way_miRBase_v21_hg38Tohg19_f1_f23.txt', header=F, row.names=1)
colnames(data) <- c('1','2','3','4','5','6','7','8','9','10',
                    '11','12','13','14','15','16','17','18','19',
                    '20','21','22')
df <- melt(data)

p.boxplot <- ggplot(df,aes(x=variable,y=value))
p.boxplot <- p.boxplot + geom_boxplot()
p.boxplot <- p.boxplot + ggtitle("miRNA conservation score[1-22nt] - phyloP")
p.boxplot <- p.boxplot + xlab("nt") + ylab("Conservation score(phyloP)")
p.boxplot <- p.boxplot + theme(plot.title = element_text(size=25),
                               axis.text.x = element_text(size=15),
                               axis.text.y = element_text(size=15),
                               axis.title.x = element_text(size=20),
                               axis.title.y = element_text(size=20))
plot(p1)
