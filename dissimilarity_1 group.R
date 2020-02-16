# Dissimilarity test for one group 

rm(list=ls())
#setwd("E:/桌面/Go-PCR实验/5.8 哀牢山森林土样本/数据分析6.5/做了ITSx/uparse")
setwd("E:/桌面/test.data/")
suppressMessages(library(vegan))

fg.dat = read.table("Galaxy520-[resample_Galaxy292-normalized_uparse_otu_ALS_comm2.txt].txt", sep="\t", row.names=1, header=T)
list.dat = read.table("group_1.txt", sep="\t", row.names=1 )

# fg.dat<-read.table(file="Galaxy73-[UPARSE_otu_table.txt].tabular",sep="\t",header=T,row.names=1)
# list.dat = read.table("first+second+groupmaking.txt", sep="\t", row.names=1 )

first.group = list.dat[, 1]
grp1 = unique(first.group);grp1
 
fg.dat2 = fg.dat[, 1:ncol(fg.dat)]#1 is for without taxonomy; 8 is for taxonomy;
fg.dat2[is.na(fg.dat2)] = 0
sample = colnames(fg.dat2)

# group cutting
samp.fg = colnames(fg.dat)
samp.env= rownames(list.dat)
my.fg = match(samp.env, samp.fg)    ##向前一个即OTU对齐。出现NA是组中的样本比OTU还多。
#my.fg2 = match(samp.fg,samp.env)  
fg.dat2 = fg.dat2[,my.fg]
rowS = rowSums(fg.dat2>0)
valid.row = which(rowS>0)
fg.dat2 = fg.dat2[valid.row, ]

#=======================================================
# for Dissimilarity test among group1
grp<-list()
for (i in c(1:length(grp1))){
  a<-as.vector(grp1[[i]])  ##组名
  grp[[i]]<-rownames(list.dat)[which(first.group==a)]
}
grp
names(grp)<-grp1 ; names(grp)
samp<-colnames(fg.dat2) ; samp
mrpp.re = matrix(0, nrow=length(grp), ncol=length(grp))
ado.re = matrix(0, nrow=length(grp), ncol=length(grp))
ano.re = matrix(0, nrow=length(grp), ncol=length(grp))

for (x in c(1:(length(grp)-1))) {
  for(y in c((x+1):length(grp))){
    list1 = grp[[x]]
    list2 = grp[[y]]
    
    # ?pmatch.pmatch(x, table) 
    col1 = pmatch(list1, samp)  ##在samp中找出list1对应的列数
    col2 = pmatch(list2, samp)
    grp.list = c(rep(as.vector(grp1[[x]]),length(list1)),rep(as.vector(grp1[[y]]), length(list2)))
    dat = fg.dat2[, c(col1, col2)]
    #====cut empty row====
    sum1 = rowSums(dat, na.rm=T) 
    valid.row = which(sum1 > 0)
    #=====================
    dat = dat[valid.row,]
    dat[is.na(dat)] = 0
    dat1 = t(dat)  #注意转置
    
    #dat.mrpp$delta, dat.mrpp$Pvalue, dat.ano$statistic, dat.ano$signif, dat.ado$aov.tab[1,4], dat.ado$aov.tab[1,6]
    #dat.dist = vegdist(dat1, method = "jaccard",binary=TRUE)  ##bray
    dat.dist = vegdist(dat1, method = "bray")  ##bray
    ?mrpp
    dat.mrpp = mrpp(dat.dist, grp.list)  
    #str(dat.mrpp)
    mrpp.re[x, y] = dat.mrpp$Pvalue  
    mrpp.re[y, x] = dat.mrpp$delta   #上三角是P，下三角是特征值
    
    ?anosim
    dat.ano = anosim(dat.dist, grp.list) #bray
    ano.re[x, y] = dat.ano$signif  
    ano.re[y, x] = dat.ano$statistic
    
    ?adonis
    grp.vector = list(V1 = grp.list)
    #dat.ado = adonis(dat1 ~ V1, data=grp.vector, method = "jaccard",binary=TRUE)
    dat.ado = adonis(dat1 ~ V1, data=grp.vector, method = "bray")
    ado.re[y, x] = dat.ado$aov.tab[1,4]
    ado.re[x, y] = dat.ado$aov.tab[1,6] 
    print(dat.ado)
  }
}
colnames(mrpp.re) = rownames(mrpp.re) <- grp1
mrpp.re
colnames(ado.re) = rownames(ado.re) <- grp1
ado.re
colnames(ano.re) = rownames(ano.re) <- grp1
ano.re
write.table(mrpp.re,file="mrpp.txt",sep="\t",col.names=NA)
write.table(ado.re,file="adonis.txt",sep="\t",col.names=NA)
write.table(ano.re,file="anosim.txt",sep="\t",col.names=NA)

#===============================================================
# do dissimilarity for whole dataset based on the grp1 grouping profile
fg.dat2<-fg.dat2[,as.vector(unlist(grp))]
fg.dat2[is.na(fg.dat2)] = 0
grp.list = c()
for(i in 1:length(grp)){
  grp.list = c(grp.list, rep(paste("grp",i,sep=""),length(grp[[i]])))
}
grp.list

report=c()
grp.vector = list(V1 = grp.list)

dat.dist = vegdist(t(fg.dat2), method = "jaccard", binary=TRUE) #"bray"

dat.mrpp = mrpp(dat.dist, grp.list) 

dat.ano = anosim(dat.dist, grp.list) 

dat.ado = adonis(t(fg.dat2) ~ V1, data=grp.vector, method = "jaccard", binary=TRUE)
#dat.ado = adonis(t(fg.dat2) ~ V1, data=grp.vector, method = "bray")

report = cbind(report, c(dat.mrpp$delta, dat.mrpp$Pvalue, dat.ano$statistic, dat.ano$signif, dat.ado$aov.tab[1,4], dat.ado$aov.tab[1,6]))
rownames(report) <- c("MRPP.delta","MRPP.P","ANOSIM.r","ANOSIM.P","PERMANOVA.F","PERMANOVA.P")
colnames(report) = "Whole dataset"
report
write.table(report,file="whole_group.txt",sep="\t",col.names = NA)
