#Title: phastcons_score_heatmap
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-06-10

#library
require(ggplot2)
require(data.table)
require(reshape2)

#file_path
setwd("C:/Users/Naoto/Documents/github/MIRAGE/data/PhastCons46Ways/")
input_file_path = "C:/Users/Naoto/Documents/github/MIRAGE/data/PhastCons46Ways/phastCons46way_miRBase_v21_hg38Tohg19.txt"
output_file_path = "C:/Users/Naoto/Documents/github/MIRAGE/data/PhastCons46Ways/phastCons46way_miRBase_v21_hg38Tohg19_heatmap.png"

#read_files
label <- c('name','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28')
setnames(input_data <- fread(input_file_path,header=F)[1:1000],label)
input_data.melt <- melt(input_data)

p1 <- ggplot(input_data.melt, aes(variable,name)) 
p1 <- p1 + geom_tile(aes(fill=value),colour="white") 
p1 <- p1 + scale_fill_gradient(low="black",high="steelblue")
plot(p1)