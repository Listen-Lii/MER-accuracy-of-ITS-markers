rm(list=ls(all=TRUE))

##ggplot-diversity
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


#########################################################################PD_PCoA
####### PD

setwd("E:/桌面/test.data/")
#setwd("E:/桌面/Go-PCR实验/5.8 哀牢山森林土样本/数据分析6.5/做了ITSx/YN")
library("picante");library("ape");library("vegan")
p_load(picante,ape,vegan)
# phytree = read.tree("Galaxy33-[FastTree.nwk].nwk");phytree  ##FastTree结果文件
# summary(phytree)
# otu = read.table("Galaxy20-[resample1400_normalized_128-uparse_otu.txt].txt", header=T, row.names=1, sep="\t")
# group = read.table("group.txt", sep="\t", row.names=1 )

phytree = read.tree("Galaxy429-[FastTree.nwk].nwk");phytree  ##FastTree结果文件
otu = read.table(file="Comnine-2fold-Galaxy296-resample_uparse_otu_ALS.txt",sep="\t",header=T,row.names=1)
group = read.table("group2.txt", sep="\t", row.names=1 )

#构建有根树=====
is.rooted(phytree) ;?root # whether the tree is rooted. Tree from FastTree is unrooted.FLASE
phytree2 = root(phytree, 1, r=TRUE) ;is.rooted(phytree2) # pick the the farthest OTU as a root.
#resolve.root:	 a logical specifying whether to resolve the new root as a bifurcating node.

##修剪树和OTU===
#?prune.sample #对树进行修剪，仅包含otu中出现的物种。Prune a phylogenetic tree to include only species present in a otuunity data set or with non-missing trait data
phy.tree = prune.sample(t(otu), phytree2) #delete useless OTUs in the tree

#?match.phylo.otu #pruning and sorting the two kinds of data to match one another for subsequent analysis.
match.otu <- match.phylo.comm(phy.tree,t(otu)) # match the OTUs in the tree with the OTUs in the OTU table

#otu = c(1,2,3,4,5,6,7); tree = c(3,4,5,6,7,8,9)

str(match.otu)
otu = match.otu$comm;otu  ##已经转至了
phy = match.otu$phy;phy

##计算PD值====
# PD = pd(otu,phy, include.root=TRUE);PD;?pd
# write.table(PD, file="PD.txt",sep="\t", row.names=TRUE, col.names=NA)

#####PCoA analysis
#install.packages("GUniFrac")
library(GUniFrac) ;?GUniFrac#用于计算Unifrac距离
unifracs <- GUniFrac(otu,phy,alpha=c(0, 0.5, 1))  ;str(unifracs)

du <- unifracs$unifracs[, , "d_UW"] # 计算Unweighted UniFrac距离
dw <- unifracs$unifracs[, , "d_1"]  # # 计算weighted UniFrac距离

?pcoa  ##ape package
PCOA_unweight <- pcoa(du, correction="none", rn=NULL) #Unweighted
str(PCOA_unweight)
eig_unweight = PCOA_unweight$values
sample_unweight =  PCOA_unweight$vectors

# sink("PCoA_unweight.txt")
# PCOA_unweight
# sink()


PCOA_weight <- pcoa(dw, correction="none", rn=NULL) #weighted
str(PCOA_unweight)
eig_weight = PCOA_weight$values
sample_weight =  PCOA_weight$vectors

# sink("PCoA_weight.txt")
# PCOA_unweight
# sink()

#####cmdscale
# ?cmdscale
# pcoa_unweight = cmdscale(du, k=2, eig=T) # k is dimension, ; eig is eigenvalues
# #str(pcoa_unweight)
# #points = pcoa_unweight$points
# pcoa_weight = cmdscale(dw, k=2, eig=T) # k is dimension, ; eig is eigenvalues
# 
# 
# ##plot
# par(mfrow=c(2,1))
# x = as.matrix(sample_weight[,1:2])
# row=row.names(x)
# pc1=round(eig_weight[1,2],2);pc1
# pc2=round(eig_weight[2,2],2);pc2
# plot(x,main="weighted PCoA",xlab=paste("pc1=",pc1),ylab=paste("pc2=",pc2));text(x,row)
# points(x, pch=21, col="red", bg="red", cex=0.5)
# 
# y = as.matrix(sample_unweight[,1:2])
# row=row.names(y)
# pc1=round(eig_unweight[1,2],2);pc1
# pc2=round(eig_unweight[2,2],2);pc2
# plot(y,main="unweighted PCoA",xlab=paste("pc1=",pc1),ylab=paste("pc2=",pc2));text(y,row)
# points(y, pch=21, col="red", bg="red", cex=0.5)

#####ggplot
wei = sample_weight[,1:2]
PCoA1=round(eig_weight[1,2],2);PCoA1
PCoA2=round(eig_weight[2,2],2);PCoA2

unwei = sample_unweight[,1:2]
PCoA1=round(eig_unweight[1,2],2);PCoA1
PCoA2=round(eig_unweight[2,2],2);PCoA2
##match group and analysis results
rowgroup = rownames(group)
rowsam = rownames(unwei)  ##wei替换为unwei即可
mat = match(rowgroup,rowsam)
unwei = unwei[mat,];unwei

colnames(group) = c("Primer","Treatment","ellipse")
data.plot = cbind(unwei,group)
data.plot
#ggplot2
library(ggplot2)
rowname = row.names(unwei)

# 论文
data.plot = subset(data.plot,Primer != "BC")

p = ggplot(data.plot,aes(Axis.1,Axis.2))
p = p + geom_point(aes(colour = Primer,shape=Treatment),alpha= 1,size = 5 );p 
p = p + xlab(paste("PCoA1=",PCoA1*100,"%",sep=""))+ylab(paste("PCoA2=",PCoA2*100,"%",sep=""));p  #+ labs(title = "PCoA analysis") ;p
#p = p + geom_text(aes(label =rowname,colour = group),position = position_dodge(0.5),vjust=-1) ;p
#p = p +  plot_theme + theme(plot.title=element_text(hjust=0.5));p
p = p +theme_bw();p
#p = p +  stat_ellipse(aes(Axis.1,Axis.2,fill=group),geom="polygon", level=0.95, alpha=0.2);p  ##椭圆
#p = p +  guides(color=guide_legend("Primer"),fill=guide_legend("Treatment"));p
p=p+theme(axis.text= element_text(size=24, color="black", family  = "serif", face= "bold", vjust=0.5, hjust=0.5))+
  theme(axis.title = element_text(size=24, color="black", family  = "serif",face= "bold", vjust=0.5, hjust=0.5))+
  theme(legend.text = element_text(colour = 'black', size = 24,  family  = "serif",face = 'bold'))+
  theme(legend.title=element_blank());p

p=p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank());p
  
p = p +scale_color_discrete(breaks = c('A','B','C','BC'));p  ##改变图例顺序

# mycolor = c("#000000","#999999","#0033FF","#9900FF","#003300","#FF3300","#003333","#00FF00")
# p = p +scale_colour_manual(values=mycolor)
# p

#show.legend 不显示这一层的图例
p = p +  stat_ellipse(aes(colour = Primer, group=ellipse),level=0.95,linetype =2 ,size=0.7 ,show.legend = FALSE);p  

#ggsave("unweight.pdf", width = 4, height = 2.5,dpi = 600); ?ggsave; getwd()

library(eoffice)
f = "E:/桌面/毕业/毕业论文图片/ITS.pptx"
topptx(p,f,append=TRUE,width = 10,height = 6)
