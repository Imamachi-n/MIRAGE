#Title: phastcons_score_heatmap
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-06-10

#library
library(pheatmap)

#file_path
setwd("C:/Users/Naoto/Documents/github/MIRAGE/data/PhyloP/")
input_file_path = "C:/Users/Naoto/Documents/github/MIRAGE/data/PhastCons46Ways/phastCons46way_miRBase_v21_hg38Tohg19_f1_f23.txt"

#output_text.km4 <- "C:/Users/Naoto/Documents/github/MIRAGE/data/PhastCons46Ways/PhastCons46Ways_miRBase_v21_hg38Tohg19_kmeans4.txt"
output_text.km5 <- "C:/Users/Naoto/Documents/github/MIRAGE/data/PhastCons46Ways/PhastCons46Ways_miRBase_v21_hg38Tohg19_kmeans5.txt"
#output_text.km6 <- "C:/Users/Naoto/Documents/github/MIRAGE/data/PhastCons46Ways/PhastCons46Ways_miRBase_v21_hg38Tohg19_kmeans6.txt"
#output_text.km7 <- "C:/Users/Naoto/Documents/github/MIRAGE/data/PhastCons46Ways/PhastCons46Ways_miRBase_v21_hg38Tohg19_kmeans7.txt"

output_file_path = "C:/Users/Naoto/Documents/github/MIRAGE/data/PhastCons46Ways/PhastCons46Ways46way_miRBase_v21_hg38Tohg19_heatmap.png"

data <- read.table(input_file_path,row.names=1,header=F)
colnames(data) <- c('1','2','3','4','5','6','7','8','9','10',
                    '11','12','13','14','15','16','17','18','19',
                    '20','21','22')

#d.km4 <- kmeans(as.matrix(data),centers=4)
d.km5 <- kmeans(as.matrix(data),centers=5)
#d.km6 <- kmeans(as.matrix(data),centers=6)
#d.km7 <- kmeans(as.matrix(data),centers=7)

#result.km4 <- cbind(d.km4$cluster,data)
result.km5 <- cbind(d.km5$cluster,data)
#result.km6 <- cbind(d.km6$cluster,data)
#result.km7 <- cbind(d.km7$cluster,data)

#o<- order(result.km4[,1]) # order the last column
#result.km4 <- result.km4[o,] # order the matrix according to the order of the last column
#pheatmap( result.km4[,2:23], cluster_rows = F, cluster_cols = F)

o<- order(result.km5[,1])
result.km5 <- result.km5[o,]
pheatmap( result.km5[,2:23], cluster_rows = F, cluster_cols = F)

#o<- order(result.km6[,1])
#result.km6 <- result.km6[o,]
#pheatmap( result.km6[,2:23], cluster_rows = F, cluster_cols = F)

#o<- order(result.km7[,1])
#result.km7 <- result.km7[o,]
#pheatmap( result.km7[,2:23], cluster_rows = F, cluster_cols = F)

#write.table(file=output_text.km4,result.km4,sep="\t",quote=F)
write.table(file=output_text.km5,result.km5,sep="\t",quote=F)
#write.table(file=output_text.km6,result.km6,sep="\t",quote=F)
#write.table(file=output_text.km7,result.km7,sep="\t",quote=F)
