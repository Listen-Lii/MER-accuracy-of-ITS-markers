###正态性检验
#2019.5.10

library(hilldiv)


?shapiro.test()#p大于0.05，则认为不能拒绝为正态
?bartlett.test()#先确定正态分布再用这个。p大于0.05，则认为方差齐性
?fligner.test #非参检验，不依赖于任何分布。检验方差齐性



# ##figure 1 chimera ------------------------------------------------------


setwd("E:/桌面/ITS引物评价文章/2018.1.10 十二稿-Molcular ecology resources/")
otu.table = read.table("chimera.txt",sep="\t",header=T,row.names=1)
hierarchy.table = read.table("Hierarchy.table.txt",sep="\t",header=T)
# otu = otu.table[,c(1:4)]
# hierarchy = hierarchy.table[c(1:4),]
# div.test(otu,qvalue=0,hierarchy=hierarchy)

#12列。
sh.p = c()
for (i in 1:12){  # i=1 
  sh = shapiro.test(otu.table[,i])
  sh.p = cbind(sh.p,sh$p.value)
};sh.p

#taq vs kapa (1-2,3-4,...)
ba.p = c()
for (i in seq(1,11,2)){  # i=1 
  ba = fligner.test(otu.table[,c(i,i+1)])  #fligner.test;bartlett.test
  ba.p = cbind(ba.p,ba$p.value)
};ba.p

#primer(1-3-5,2-4-6,7-9-11,8-10-12)
ba.p <- vector(mode="numeric",length=12)
#Even
ba.p[1] = fligner.test(otu.table[,c(1,3)])$p.value  #bartlett.test;fligner.test
ba.p[2] = fligner.test(otu.table[,c(3,5)])$p.value
ba.p[3] = fligner.test(otu.table[,c(1,5)])$p.value
ba.p[4] = fligner.test(otu.table[,c(2,4)])$p.value
ba.p[5] = fligner.test(otu.table[,c(4,6)])$p.value
ba.p[6] = fligner.test(otu.table[,c(2,6)])$p.value
#stagged
ba.p[7] = fligner.test(otu.table[,c(7,9)])$p.value
ba.p[8] = fligner.test(otu.table[,c(9,11)])$p.value
ba.p[9] = fligner.test(otu.table[,c(7,11)])$p.value
ba.p[10] = fligner.test(otu.table[,c(8,10)])$p.value
ba.p[11] = fligner.test(otu.table[,c(10,12)])$p.value
ba.p[12] = fligner.test(otu.table[,c(8,12)])$p.value
ba.p
#Wilcoxon Rank Sum Test
?wilcox.test
wilcox.test(otu.table[,1],otu.table[,3],paired = F)$p.value

##
library(mvnormtest);?mshapiro.test() ##这个包只有这一个函数，进行多变量的正态性检验
mshapiro.test(t(otu))

#https://bbs.pinggu.org/thread-417275-1-1.html
#Shapiro-Wilk检验只适用于小样本场合（3≤n≤50）



# ### figure s4 soil chimera ----------------------------------------------


setwd("E:/桌面/ITS引物评价文章/2018.1.10 十二稿-Molcular ecology resources/")
otu.table = read.table("soil_chimera.txt",sep="\t",header=T,row.names=1)

#6列。
sh.p = c()
for (i in 1:6){  # i=1 
  sh = shapiro.test(otu.table[,i])
  sh.p = cbind(sh.p,sh$p.value)
};sh.p

#taq vs kapa (1-2,3-4,...)
ba.p = c()
for (i in seq(1,5,2)){  # i=1 
  ba = fligner.test(otu.table[,c(i,i+1)])  #fligner.test;bartlett.test
  ba.p = cbind(ba.p,ba$p.value)
};ba.p

#primer(1-3-5,2-4-6)
ba.p <- vector(mode="numeric",length=6)
#primer
ba.p[1] = fligner.test(otu.table[,c(1,3)])$p.value  #bartlett.test;fligner.test
ba.p[2] = fligner.test(otu.table[,c(3,5)])$p.value
ba.p[3] = fligner.test(otu.table[,c(1,5)])$p.value
ba.p[4] = fligner.test(otu.table[,c(2,4)])$p.value
ba.p[5] = fligner.test(otu.table[,c(4,6)])$p.value
ba.p[6] = fligner.test(otu.table[,c(2,6)])$p.value
ba.p



# ###figure s5 soil alpha diverisity --------------------------------------

setwd("E:/桌面/ITS引物评价文章/2018.1.10 十二稿-Molcular ecology resources/")
otu.table = read.table("alpha.diverisity.txt",sep="\t",header=T,row.names=1)

############richness
otu = otu.table[1:5,]
otu = otu[,-3]  #这一列只有两个，做shapiro的时候去掉。注意做ba的时候 还用上面一行的otu
#6列。
sh.p = c()
for (i in 1:5){  # i=1 
  sh = shapiro.test(otu[,i])
  sh.p = cbind(sh.p,sh$p.value)
};sh.p

#taq vs kapa (1-2,3-4,...)
ba.p = c()
for (i in seq(1,5,2)){  # i=1 
  ba = fligner.test(otu[,c(i,i+1)])  #fligner.test;bartlett.test
  ba.p = cbind(ba.p,ba$p.value)
};ba.p

#primer(1-3-5,2-4-6)
ba.p <- vector(mode="numeric",length=6)
#primer
ba.p[1] = fligner.test(otu[,c(1,3)])$p.value  #bartlett.test;fligner.test
ba.p[2] = fligner.test(otu[,c(3,5)])$p.value
ba.p[3] = fligner.test(otu[,c(1,5)])$p.value
ba.p[4] = fligner.test(otu[,c(2,4)])$p.value
ba.p[5] = fligner.test(otu[,c(4,6)])$p.value
ba.p[6] = fligner.test(otu[,c(2,6)])$p.value
ba.p

#######shannon
otu = otu.table[6:10,]
otu = otu[,-3]  #这一列只有两个，做shapiro的时候去掉。注意做ba的时候 还用上面一行的otu

#6列。
sh.p = c()
for (i in 1:5){  # i=1 
  sh = shapiro.test(otu[,i])
  sh.p = cbind(sh.p,sh$p.value)
};sh.p

#taq vs kapa (1-2,3-4,...)
ba.p = c()
for (i in seq(1,5,2)){  # i=1 
  ba = fligner.test(otu[,c(i,i+1)])  #fligner.test;bartlett.test
  ba.p = cbind(ba.p,ba$p.value)
};ba.p

#primer(1-3-5,2-4-6)
ba.p <- vector(mode="numeric",length=6)
#primer
ba.p[1] = fligner.test(otu[,c(1,3)])$p.value  #bartlett.test;fligner.test
ba.p[2] = fligner.test(otu[,c(3,5)])$p.value
ba.p[3] = fligner.test(otu[,c(1,5)])$p.value
ba.p[4] = fligner.test(otu[,c(2,4)])$p.value
ba.p[5] = fligner.test(otu[,c(4,6)])$p.value
ba.p[6] = fligner.test(otu[,c(2,6)])$p.value
ba.p
