#########################################################################chao
rm(list=ls())
setwd("E:/桌面/test.data/")
suppressMessages(library(vegan))
##得到chao值后画图
#注意原始文件需要改名字
x<-read.table(file="Galaxy461-[chao.txt].txt",sep="\t",header=T,row.names=1)
group = read.table("group3.txt", sep="\t", row.names=1 )

library(ggplot2)

plot_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=7),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=7),
                   text=element_text(family="sans", size=7)
) + theme_bw()



sam=x

rowgroup = rownames(group)
rowsam = rownames(sam)  
mat = match(rowgroup,rowsam)  ##%in%
sam = sam[mat,]#;sam

colnames(group) = "group"
data.plot = cbind(sam,group) #data.plot

library(Rmisc)
chao_sd <- summarySE(data.plot, measurevar="sam", groupvars="group")

p = ggplot(chao_sd, aes(x=group, y=sam, fill=group)) + 
  geom_bar(position=position_dodge(), stat="identity",size=0.3) + 
  geom_errorbar(aes(ymin=sam-sd, ymax=sam+sd), width=.2, size =.3, position=position_dodge(.9)) +
  plot_theme+ theme(legend.position="None",axis.title.x = element_blank())+labs(y="Chao value")+
  theme(axis.text= element_text(size=16, color="black", family  = "serif", face= "bold", vjust=0.5, hjust=0.5))+
  theme(axis.title = element_text(size=16, color="black", family  = "serif",face= "bold", vjust=0.5, hjust=0.5))+
  ylim(c(0,3500))
p=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p


p
#ggsave(paste("chao.pdf", sep=""), p, width = 8, height = 5)

