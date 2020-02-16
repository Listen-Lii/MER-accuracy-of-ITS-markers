#########################################################################alpha diversity
rm(list=ls())
setwd("E:/桌面/test.data/")
suppressMessages(library(vegan))

library(ggplot2)
plot_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=14),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=7),
                   text=element_text(family="sans", size=7)
) +theme_bw()



x<-read.table(file="Comnine-2fold-Galaxy296-resample_uparse_otu_ALS.txt",sep="\t",header=T,row.names=1)
x[is.na(x)] = 0
x = t(x)

Shannon <- diversity(x)
Inv_Simpson <- diversity(x, "inv")
S <- specnumber(x)
Pielou_evenness <- Shannon/log(S)
Simpson_evenness <- Inv_Simpson/S
Observed_richness <- colSums(t(x)>0)
report1 = cbind(Shannon, Inv_Simpson,Observed_richness, Pielou_evenness, Simpson_evenness)

#write.table(report1, "com-alpha.txt", sep="\t", col.names = NA)

##########ggplot2 barplot
library(ggplot2)
group = read.table("group3.txt", sep="\t", row.names=1 )
sam=report1

rowgroup = rownames(group)
rowsam = rownames(sam)        ##deblur注意改样本名称
mat = match(rowgroup,rowsam)  ##%in%
sam = sam[mat,]#;sam

colnames(group) = "group"
data.plot = cbind(sam,group) #data.plot

library(Rmisc)
shannon_sd <- summarySE(data.plot, measurevar="Shannon", groupvars="group")
Inv_Simpson_sd <- summarySE(data.plot, measurevar="Inv_Simpson", groupvars="group")
Observed_richness_sd <- summarySE(data.plot, measurevar="Observed_richness", groupvars="group")
Pielou_evenness_sd <- summarySE(data.plot, measurevar="Pielou_evenness", groupvars="group")

##散点图
# ggplot(shannon_SE, aes(x=group, y=Shannon, colour=group)) + 
#   geom_errorbar(aes(ymin=Shannon-sd, ymax=Shannon+sd), width=.1) + 
#   geom_line() + geom_point(size = 3) 

##柱状图
# windowsFonts()
# windowsFonts(HEL=windowsFont("Helvetica CE 55 Roman"),
#              RMN=windowsFont("Times New Roman"),
#              ARL=windowsFont("Arial"))
# windowsFonts(myFont1=windowsFont("Times New Roman"),myFont2=windowsFont("华文行楷"))

p1 = ggplot(shannon_sd, aes(x=group, y=Shannon, fill=group)) + 
  geom_bar(position=position_dodge(), stat="identity",size=0.3) + 
  geom_errorbar(aes(ymin=Shannon-sd, ymax=Shannon+sd), width=.2, size =.3, position=position_dodge(.9)) +
  plot_theme+ theme(legend.position="None",axis.title.x = element_blank())+
  theme(axis.text= element_text(size=16, color="black", family  = "serif", face= "bold", vjust=0.5, hjust=0.5))+
  theme(axis.title = element_text(size=16, color="black", family  = "serif",face= "bold", vjust=0.5, hjust=0.5))+
  ylim(c(0,7.3))
p1=p1+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p1



library(ggpubr)
my_comparisons1 <- list(c("gITS7/ITS4R-e", "gITS7/ITS4R-c"), c("5.8S-Fun/ITS4-Fun-c", "5.8S-Fun/ITS4-Fun-e"), c("ITS1F/ITS2-c", "ITS1F/ITS2-e"))
                        
my_comparisons2 <- list(c("ITS1F/ITS2-c", "5.8S-Fun/ITS4-Fun-c" ),c("ITS1F/ITS2-c","gITS7/ITS4R-c"),c("gITS7/ITS4R-c", "5.8S-Fun/ITS4-Fun-c" ))

my_comparisons3 <- list(c("gITS7/ITS4R-e", "ITS1F/ITS2-e"),c("gITS7/ITS4R-e", "5.8S-Fun/ITS4-Fun-e"),c("5.8S-Fun/ITS4-Fun-e", "ITS1F/ITS2-e"))


#compare_means(Shannon~group,shannon_sd,method = "kruskal.test")
p1+ stat_compare_means(comparisons=my_comparisons1,label.y = c(6, 6, 6))+ 
  stat_compare_means(comparisons=my_comparisons2,label.y = c(6.9, 7, 7.1))+
  stat_compare_means(comparisons=my_comparisons3,label.y = c(6.4, 6.5, 6.6))+

  stat_compare_means(label = "p.signif", ref.group = ".all.", hide.ns = TRUE) 


p2 = ggplot(Inv_Simpson_sd, aes(x=group, y=Inv_Simpson, fill=group)) + 
  geom_bar(position=position_dodge(), stat="identity",size=0.3) + 
  geom_errorbar(aes(ymin=Inv_Simpson-sd, ymax=Inv_Simpson+sd), width=.2, size =.3, position=position_dodge(.9)) +
  plot_theme+ theme(legend.position="None",axis.title.x = element_blank())+
  theme(axis.text= element_text(size=16, color="black", family  = "serif", face= "bold", vjust=0.5, hjust=0.5))+
  theme(axis.title = element_text(size=16, color="black", family  = "serif",face= "bold", vjust=0.5, hjust=0.5))+
  ylim(c(0,55))
p2=p2+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

p2
p2+ stat_compare_means(comparisons=my_comparisons)

  
p3 = ggplot(Observed_richness_sd, aes(x=group, y=Observed_richness, fill=group)) + 
  geom_bar(position=position_dodge(), stat="identity",size=0.3) + 
  geom_errorbar(aes(ymin=Observed_richness-sd, ymax=Observed_richness+sd), width=.2, size =.3, position=position_dodge(.9)) +
  plot_theme+ theme(legend.position="None",axis.title.x = element_blank())+
theme(axis.text= element_text(size=16, color="black", family  = "serif", face= "bold", vjust=0.5, hjust=0.5))+
  theme(axis.title = element_text(size=16, color="black", family  = "serif",face= "bold", vjust=0.5, hjust=0.5))+
  ylim(c(0,2300))
p3=p3+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

p3

p3+ stat_compare_means(comparisons=my_comparisons)

p4 = ggplot(Pielou_evenness_sd, aes(x=group, y=Pielou_evenness, fill=group)) + 
  geom_bar(position=position_dodge(), stat="identity",size=0.3) + 
  geom_errorbar(aes(ymin=Pielou_evenness-sd, ymax=Pielou_evenness+sd), width=.2, size =.3, position=position_dodge(.9)) +
  plot_theme+ theme_bw()+ theme(legend.position="None",axis.title.x = element_blank())
p4

#1.分面
library(gridExtra)
grid.arrange(p1, p2, p3, p4)