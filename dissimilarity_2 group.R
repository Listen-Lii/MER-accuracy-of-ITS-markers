# Dissimilarity test for two or more groups

rm(list=ls())
setwd("E:/桌面/test.data/")
suppressMessages(library(vegan))

x<-read.table(file="Comnine-2fold-Galaxy296-resample_uparse_otu_ALS.txt",sep="\t",header=T,row.names=1)
list.dat = read.table("group_7F_adonis.txt", sep="\t", row.names=1 )
fg.dat = x

first.group = list.dat[, 1]
grp1 = unique(first.group) ;grp1
second.group = list.dat[, 2]
grp2 = unique(second.group) ;grp2

 
fg.dat2 = fg.dat[, 1:ncol(fg.dat)]#1 is for without taxonomy; 8 is for taxonomy;
sample = colnames(fg.dat2)

# match
samp.fg = colnames(fg.dat)
samp.env= rownames(list.dat)
my.fg = match(samp.env, samp.fg)  # ?match match(x, table)
fg.dat2 = fg.dat2[,my.fg]
rowS = rowSums(fg.dat2>0)
valid.row = which(rowS>0)
fg.dat2 = fg.dat2[valid.row, ]

# for adonis,dissimilarity between group1 and 2
  list.dat.inter = list(V1=list.dat[,1],V2=list.dat[,2])
#  dat.ado = adonis(t(fg.dat2) ~ V1*V2, data=list.dat.inter, method = "jaccard", binary = T) 
 dat.ado = adonis(t(fg.dat2) ~ V1*V2, data=list.dat.inter, method = "bray")
 
 list.dat.inter = list(V1=list.dat[,1])
 dat.ado = adonis(t(fg.dat2) ~ V1, data=list.dat.inter, method = "bray")
# y ~model 是一种特定的格式， 表示以 y 为响应变量， 模型为 model。 model 中的变量由 + 来连接， 或者由: 来表示变量间的 \交互作用"。a*b=a+b+a:b
#  (a+b+c)^2表示(a+b+c)*(a+b+c)。
# -表示去掉。(a + b + c)^2 - a : b 表示 a + b + c + b : c + a : c
  
  dat.ado
  str(dat.ado)
  dat.ado$aov.tab
  
  report=c()
  report1 = cbind(dat.ado$aov.tab[1,4],dat.ado$aov.tab[1,5],dat.ado$aov.tab[1,6])#F,R,p
  report2 = cbind(dat.ado$aov.tab[2,4],dat.ado$aov.tab[2,5],dat.ado$aov.tab[2,6])
  report3 = cbind(dat.ado$aov.tab[3,4],dat.ado$aov.tab[3,5],dat.ado$aov.tab[3,6])
  report=rbind(report1,report2,report3)
  colnames(report)<-c("F.model","R2","P-value")
  rownames(report)<-c("Group1","Group2","Group1*Group2")
  report
  
  write.table(report,file="adonis_with_comm_jaccard.txt",sep="\t",col.names = NA)
  
  
# for each grp1, testing dissimilarity of group2
  report = c()
  for(x in 1:length(grp1)){
    first.col = which(first.group == grp1[x])
    grp.list = second.group[first.col]
    dat = fg.dat2[,first.col]  
    # delete empty rows
    rsum1 = rowSums(dat)
    tempCK1 = which(rsum1==0)
    if(length(tempCK1)!=0) {dat = dat[-tempCK1,]}
    # calculate dissimilarity
    dat1 = t(dat)
    dat.dist = vegdist(dat1, method = "jaccard", binary=T)
   #dat.dist = vegdist(dat1, method = "bray")
    
    dat.mrpp = mrpp(dat.dist, grp.list)
    
    dat.ano = anosim(dat.dist, grp.list)

    grp.vector = list(V1 = grp.list)
    dat.ado = adonis(dat1 ~ V1, data=grp.vector, method = "jaccard", binary=T)   
  # dat.ado = adonis(dat1 ~ V1, data=grp.vector, method = "bray")  
    report = rbind(report, c(dat.mrpp$delta, dat.mrpp$Pvalue, dat.ano$statistic, dat.ano$signif, dat.ado$aov.tab[1,4], dat.ado$aov.tab[1,6]))
  }
  ##如果数据集太小，会出现'nperm' >= set of all permutations: complete enumeration.Set of permutations < 'minperm'. Generating entire set.
  ##requested more permutations (default: permutations=999) than the number of theoretically possible permutations with your data set, so anosim makes a complete enumeration, i.e. checks all possible combinations.
  
  colnames(report) = c("MRPP.delta","MRPP.P","ANOSIM.r","ANOSIM.P","PERMANOVA.F","PERMANOVA.P")
  rownames(report) = grp1
  report

  write.table(report,file="dis_within_group1.txt",sep="\t",col.names=NA)
###	

# for each grp2, testing dissimilarity of group1
#交换group1和2
  first.group = list.dat[, 2]
  grp1 = unique(first.group) ;grp1
  second.group = list.dat[, 1]
  grp2 = unique(second.group) ;grp2 
#重新从51行开始运行。
  write.table(report,file="dis_within_group2.txt",sep="\t",col.names=NA)
  