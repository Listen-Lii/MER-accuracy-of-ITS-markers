##2019.5.22
#环境样本门水平累积柱状图
getwd()
setwd("E:/桌面/ITS引物评价文章/2018.1.10 十二稿-Molcular ecology resources/文章修改文件/phylum/")
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
plot_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=20),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=20),
                   text=element_text(family="sans", size=20)
)
##soil.txt and endo.txt
x<-read.table(file="endo.txt",sep="\t",header=T)  
x
#长宽转化
#data = gather(x,Phylum,percentage,Ascomycota:Unclassified)
data = gather(x,Phylum,percentage,Ascomycota:Zygomycota)
data

display.brewer.all()
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(data, aes(x=primer, y=percentage, fill=Phylum)) +
  facet_grid(. ~ treatment,scales="free") +
  theme(legend.position="right")+
  geom_bar(position="stack", stat="identity") +
  xlab("")+ylab("Percentage (%)")+
  #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #scale_fill_manual(values=cbPalette)+
  scale_fill_brewer(palette="Set3")+
  theme(axis.text= element_text(size=20, color="black", family  = "serif", face= "bold", vjust=0.5, hjust=0.5))+
  theme(axis.title = element_text(size=20, color="black", family  = "serif",face= "bold", vjust=0.5, hjust=0.5))+
  plot_theme

