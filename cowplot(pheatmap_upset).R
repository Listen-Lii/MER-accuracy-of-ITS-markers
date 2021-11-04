########################pheatmap
######参考pheatmap.R
rm(list=ls())
setwd("E:/桌面/已备份/ITS引物评价文章MER/12.2第6稿/")
genus = read.table(file="Galaxy421-(Taxonomy_summary_at_7_level_for_abundance_count).txt",sep="\t",row.names=1,header = T) 
genus = read.table(file="Galaxy421-biyelunwen.txt",sep="\t",row.names=1,header = T) 
head(genus)

# setwd("E:/桌面/已备份/ITS引物评价文章MER/2018.1.10 十二稿-Molcular ecology resources/文章修改文件/end_genus/")
# genus = read.table(file="heatmap.txt",sep="\t",row.names=1,header = T) 
# head(genus)
# genus

library(pheatmap) #packageVersion("pheatmap")  ##1.0.10
#pheatmap(genus)

p1=pheatmap(genus, scale = "row", clustering_distance_cols = "correlation",fontsize_col = 12,fontsize_row =12)  
p1
#p1 = pheatmap(genus, scale = "row", clustering_distance_rows = "correlation")  ##行为物种,按照行标准化，并按照列（样本）进行聚类

library(eoffice)
f = "E:/桌面/毕业/毕业论文图片/ITS.pptx"
topptx(p1,f,append=TRUE,width = 10,height = 6)
#ggsave(p ,filename = "pheatmap with upset.pdf")
# Fix cell sizes and save to file with correct size
# pheatmap(genus, cellwidth = 15, cellheight = 12, main = "Example heatmap")
# #输出到桌面
# pheatmap(genus, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "test.pdf")

########################upset
#####参考Upset_venn.R及Upset_ggplot.R
require(UpSetR)
mut <- read.csv("genus.csv" , header = T, sep = ",") ##注意，列是样本，行是物种。一定要先转化为0和1再做。
##注意不要数字开头，不能有-和_，_会显示为点

#upset(mut,order.by = "freq", empty.intersections = "on")
#upset(mu, order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))

p2 = upset(mut, sets = c("A_e", "A_c", "B_e", "B_c", "C_e", "C_c","BC_e","BC_c"), sets.bar.color = "red",order.by = "freq", empty.intersections = "on")
p2
require(ggplotify)
up2 <- as.ggplot(p2)

# up1 <- as.ggplot(p1)
# ?as.grob
# class(p1)  #pheatmap 不被支持转化为ggplot，只能用拼图了
# require(ggimage)
# g <- p1 + geom_subview(subview = up2 + theme_void(), x=.7, y=.7, w=.6, h=.6)
# g

##拼图.不是ggplot要加$gtable
require(cowplot)
plot_grid(p1$gtable, up2, nrow=2, labels=c("A", "B"))


mut <- read.table("upset.txt" , header = T, sep = "\t",row.names = 1)
head(mut)

library(vegan)
mut2 = decostand(mut,method="pa")
p2 = upset(mut2, sets = c("A_e", "A_c", "B_e", "B_c", "C_e", "C_c"), sets.bar.color = "red",order.by = "freq", empty.intersections = "on")

upset(mut2)
