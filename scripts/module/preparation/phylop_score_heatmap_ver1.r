#Title: phastcons_score_heatmap
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-06-10

#library
#require(ggplot2)
#require(data.table)
#require(reshape2)
library(gplots)

#file_path
setwd("C:/Users/Naoto/Documents/github/MIRAGE/data/PhyloP/")
input_file_path = "C:/Users/Naoto/Documents/github/MIRAGE/data/PhyloP/phyloP46way_miRBase_v21_hg38Tohg19_f1_f23.txt"
#input_file_path = "C:/Users/Naoto/Documents/github/MIRAGE/data/PhastCons46Ways/phastCons46way_miRBase_v21_hg38Tohg19_f1_f23.txt"
output_file_path = "C:/Users/Naoto/Documents/github/MIRAGE/data/PhyloP/phyloP46way_miRBase_v21_hg38Tohg19_heatmap.png"
#output_file_path = "C:/Users/Naoto/Documents/github/MIRAGE/data/PhastCons46Ways/PhastCons46Ways46way_miRBase_v21_hg38Tohg19_heatmap.png"

#read_files
#label <- c('name','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28')
#setnames(input_data <- fread(input_file_path,header=F)[100:200],label)
#input_data.melt <- melt(input_data)

#p1 <- ggplot(input_data.melt, aes(variable,name)) 
#p1 <- p1 + geom_tile(aes(fill=value),colour="white") 
#p1 <- p1 + scale_fill_gradient(low="white",high="steelblue")
#plot(p1)

input_file <- read.table(input_file_path,row.names=1,header=F)

heatmap(
  as.matrix(input_file),               # Data
  main = "Heat colors",     # Title
  Rowv = TRUE,              # clustering(y-axis)
  Colv = NA,                # clustering(x-axis)
  distfun = dist,           # 
  hclustfun = hclust,       # 
  col = heat.colors(256)    # color
)

heatmap.2(
  as.matrix(input_file),
  scale = "row",              # ?
  dendrogram = "none",        # 系統樹の描画を指定（both, row, column, none）
  Rowv = FALSE,               # dendrogramにboth,rowを指定した時にTRUEにする必要があります
  Colv = FALSE,               # dendrogramにboth,columnを指定した時にTRUEにする必要があります
  col = redgreen(256),
  #key = TRUE,                 # スケールを表示
  #density.info = "density",   # スケールバーに密度をグラフに示す
  main = "Density none"
)

heatmap.2(
  as.matrix(input_file), scale = "row",
  dendrogram = "both", Rowv = TRUE, Colv = TRUE,
  trace = "none",
  col = redgreen(256),
  main = "Density none"       # グラフにヒストグラムを表示させない
)