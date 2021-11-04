library(UpSetR)

##物种venn图
setwd("E:/桌面/test.data")
mut <- read.csv("genus.csv" , header = T, sep = ",") ##注意，列是样本，行是物种。一定要先转化为0和1再做。
##注意不要数字开头，不能有-和_，_会显示为点
head(mut)
upset(mut,order.by = "freq", empty.intersections = "on")

#upset(mu, order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))

upset(mut, sets = c("A_e", "A_c", "B_e", "B_c", "C_e", "C_c","BC_e","BC_c"), sets.bar.color = "red",order.by = "freq", empty.intersections = "on")

###2019.1.8
#一些样本中共有的线变成红色
##mb.ratio：控制上方条形图以及下方点图的比例。#order.by：如何排序，这里 freq 表示从大到小排序展示，其他选项有 degree 以及先按 freq 再按 degree 排序。
p1 = upset(mut, sets = c("A_e", "A_c", "B_e", "B_c", "C_e", "C_c","BC_e","BC_c"),
           # number.angles = 30, 
           point.size = 2, line.size = 1,
           mainbar.y.label = "OTU", sets.x.label = "OTU Per Treatment",
           text.scale = c(2, 2, 1.5,1.5, 1.5, 1.5),mb.ratio = c(0.7, 0.3),
           order.by = "freq",keep.order = TRUE,
           queries = list(list(query = intersects, params = 
                                 list("A_e", "A_c", "B_e", "B_c", "C_e", "C_c","BC_e","BC_c"), color = "red", active = T)))

p1


### 毕业论文
upset(mut, sets = c("A_e", "A_c", "B_e", "B_c", "C_e", "C_c"), sets.bar.color = "red",order.by = "freq", empty.intersections = "on")

###2019.1.8
#一些样本中共有的线变成红色
##mb.ratio：控制上方条形图以及下方点图的比例。#order.by：如何排序，这里 freq 表示从大到小排序展示，其他选项有 degree 以及先按 freq 再按 degree 排序。
p1 = upset(mut, sets = c("A_e", "A_c", "B_e", "B_c", "C_e", "C_c"),
           # number.angles = 30, 
           point.size = 2, line.size = 1,
           mainbar.y.label = "OTU", sets.x.label = "OTU Per Treatment",
           text.scale = c(2, 2, 1.5,1.5, 1.5, 1.5),mb.ratio = c(0.7, 0.3),
           order.by = "freq",keep.order = TRUE,
           queries = list(list(query = intersects, params = 
                                 list("A_e", "A_c", "B_e", "B_c", "C_e", "C_c"), color = "red", active = T)))

p1
library(eoffice)
f = "E:/桌面/毕业/毕业论文图片/ITS.pptx"
topptx(p1,f,append=TRUE,width = 8,height = 6)
