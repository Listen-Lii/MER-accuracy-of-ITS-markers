##2019.5.22
##mismatch 累积柱状图
##2019.5.29 调整纵坐标mismatch顺序。
setwd("E:/桌面/ITS引物评价文章/2018.1.10 十二稿-Molcular ecology resources/文章修改文件/mismatch/")
x<-read.table(file="mis.txt",sep="\t",header=T)  
x

library(tidyverse)
E_A = x %>% filter(Type == "E",primer=="A") %>%
  group_by(mismatch)%>%
  summarize(mean = mean(percentage),sd = sd(percentage)) 
E_A = cbind(Type=rep("E_A",3),E_A,for_sd=unlist(c(100,E_A[2,2],E_A[2,2]+E_A[3,2])))

E_B = x %>% filter(Type == "E",primer=="B") %>%
  group_by(mismatch) %>% 
  summarize(mean = mean(percentage),sd = sd(percentage)) 
E_B = cbind(Type=rep("E_B",3),E_B,for_sd=unlist(c(100,E_B[2,2],E_B[2,2]+E_B[3,2])))

E_C = x %>% filter(Type == "E",primer=="C") %>%
  group_by(mismatch) %>% 
  summarize(mean = mean(percentage),sd = sd(percentage)) 
E_C = cbind(Type=rep("E_C",3),E_C,for_sd=unlist(c(100,E_C[2,2],E_C[2,2]+E_C[3,2])) ) 

EK_A = x %>% filter(Type == "EK",primer=="A") %>%
  group_by(mismatch) %>% 
  summarize(mean = mean(percentage),sd = sd(percentage)) 
EK_A = cbind(Type=rep("EK_A",3),EK_A,for_sd=unlist(c(100,EK_A[2,2],EK_A[2,2]+EK_A[3,2]))  ) 

EK_B = x %>% filter(Type == "EK",primer=="B") %>%
  group_by(mismatch) %>% 
  summarize(mean = mean(percentage),sd = sd(percentage)) 
EK_B = cbind(Type=rep("EK_B",3),EK_B,for_sd=unlist(c(100,EK_B[2,2],EK_B[2,2]+EK_B[3,2])) )

EK_C = x %>% filter(Type == "EK",primer=="C") %>%
  group_by(mismatch) %>% 
  summarize(mean = mean(percentage),sd = sd(percentage)) 
EK_C = cbind(Type=rep("EK_C",3),EK_C,for_sd=unlist(c(100,EK_C[2,2],EK_C[2,2]+EK_C[3,2])) )  

S_A = x %>% filter(Type == "S",primer=="A") %>%
  group_by(mismatch) %>% 
  summarize(mean = mean(percentage),sd = sd(percentage)) 
S_A = cbind(Type=rep("S_A",3),S_A,for_sd=unlist(c(100,S_A[2,2],S_A[2,2]+S_A[3,2]))  )

S_B = x %>% filter(Type == "S",primer=="B") %>%
  group_by(mismatch) %>% 
  summarize(mean = mean(percentage),sd = sd(percentage)) 
S_B = cbind(Type=rep("S_B",3),S_B,for_sd=unlist(c(100,S_B[2,2],S_B[2,2]+S_B[3,2]))   )

S_C = x %>% filter(Type == "S",primer=="C") %>%
  group_by(mismatch) %>% 
  summarize(mean = mean(percentage),sd = sd(percentage)) 
S_C = cbind(Type=rep("S_C",3),S_C,for_sd=unlist(c(100,S_C[2,2],S_C[2,2]+S_C[3,2]))  ) 

SK_A = x %>% filter(Type == "SK",primer=="A") %>%
  group_by(mismatch) %>% 
  summarize(mean = mean(percentage),sd = sd(percentage)) 
SK_A = cbind(Type=rep("SK_A",3),SK_A,for_sd=unlist(c(100,SK_A[2,2],SK_A[2,2]+SK_A[3,2]))  ) 

SK_B = x %>% filter(Type == "SK",primer=="B") %>%
  group_by(mismatch) %>% 
  summarize(mean = mean(percentage),sd = sd(percentage)) 
SK_B = cbind(Type=rep("SK_B",3),SK_B,for_sd=unlist(c(100,SK_B[2,2],SK_B[2,2]+SK_B[3,2]))  )  

SK_C = x %>% filter(Type == "SK",primer=="C") %>%
  group_by(mismatch) %>% 
  summarize(mean = mean(percentage),sd = sd(percentage)) 
SK_C = cbind(Type=rep("SK_C",3),SK_C,for_sd=unlist(c(100,SK_C[2,2],SK_C[2,2]+SK_C[3,2])) )    

##################
da.ta = rbind(E_A,E_B,E_C,EK_A,EK_B,EK_C,S_A,S_B,S_C,SK_A,SK_B,SK_C)
data = cbind(facet=c(rep("Even community",18),rep("Staggered community",18)),da.ta)
data
data$level =rep(c("A","C","B"),12)

data2 = data %>% arrange(Type,level);data2

library(ggplot2)
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
ggplot(data, aes(x=Type, y=mean, fill=level)) +
  facet_grid(. ~ facet,scales="free") +
  theme(legend.position="right")+
  geom_bar(position="stack", stat="identity") +
  geom_errorbar(aes(ymin=for_sd+sd, ymax=for_sd-sd),
                width=.1,                    
                position="identity") +
  xlab("")+ylab("Percentage (%)")+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #labs(fill="Mismatch") + scale_fill_discrete(breaks=c(">1","1","0"))+
  theme(axis.text= element_text(size=20, color="black", family  = "serif", face= "bold", vjust=0.5, hjust=0.5))+
  theme(axis.title = element_text(size=20, color="black", family  = "serif",face= "bold", vjust=0.5, hjust=0.5))+
  plot_theme


########显著性矩阵
setwd("E:/桌面/ITS引物评价文章/2018.1.10 十二稿-Molcular ecology resources/mismatch/")
x<-read.table(file="t.txt",sep="\t",header=T)  
x
?t.test
p = c()
for (i in 1:ncol(x)){ #i=2
  p1 = c()
  for (j in 1:ncol(x)){
    t = t.test(x[,i],x[,j],paired=F)
    p1 = rbind(p1,t$p.value)
  }
  p = cbind(p,p1)
}
p
rownames(p)=colnames(x)
colnames(p)=rownames(p)
p

#
library(corrplot);?corrplot
corrplot(p,method="color",type="lower",is.corr=F,diag=F,
         p.mat = p,insig = "label_sig",col="#D1E5F0",outline=T,#矩形框
         tl.col = "black",#文本颜色
         
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "red")


# for (i in 1:(ncol(x)-1)){ #i=2   ##每次循环少一个，最后不是矩阵形式。得用上面的出矩阵。
#   p1 = c()
#   for (j in (i+1):ncol(x)){
#     t = t.test(x[,i],x[,j],paired=F)
#     p1 = rbind(p1,t$p.value)
#   }
#   p = cbind(p,p1)
# }
# p


###########资料
DF <- data.frame(sex=c("men","women","men","women","men","women"),
                 proportion=c(0.33,0.32,0.24,0.29,0.12,0.16),
                 ci_l=c(0.325,0.322,0.230,0.284,0.114,0.155),
                 ci_u=c(0.339,0.316,0.252,0.311,0.130,0.176),
                 year=c(2008,2008,2013,2013,2013,2013),
                 response=c("Yes","Yes","Yes, entire the journey","Yes, entire the journey","Yes, part of the journey","Yes, part of the journey")
)
ggplot(DF, aes(x=factor(year), y=proportion, fill=response)) +
  facet_grid(. ~ sex) +
  theme(legend.position="none")+
  geom_bar(position="stack", stat="identity") +
  geom_errorbar(aes(ymin=ci_l, ymax=ci_u),
                width=.2,                    
                position="identity") 

DF$ci_l[DF$response == "Yes, part of the journey"] <- with(DF,ci_l[response == "Yes, part of the journey"] +
                                                             ci_l[response == "Yes, entire the journey"])

DF$ci_u[DF$response == "Yes, part of the journey"] <- with(DF,ci_u[response == "Yes, part of the journey"] +
                                                             ci_u[response == "Yes, entire the journey"])

ggplot(DF, aes(x=factor(year), y=proportion)) +
  facet_grid(. ~ sex) +
  geom_bar(stat="identity",aes(fill=response)) +
  geom_errorbar(aes(ymin= ci_l, 
                    ymax= ci_u),
                width=.2,                    # Width of the error bars
                position="identity")                                                             

DF$vadj <- c(rep(0,2), rep(c(0,1,0), each=2) * DF$proportion)[1:6]

ggplot(DF, aes(x=factor(year), y=proportion, fill=response)) +
  facet_grid(. ~ sex) + geom_bar(stat='identity') + 
  geom_errorbar( aes(ymin=ci_l+vadj, ymax=ci_u+vadj), width=.2)


library(RCurl)
table.text <- getURL("https://raw.githubusercontent.com/richardbeare/ggplotissues/master/barplot_problem.csv")
dat <- read.csv(text=table.text)
#dat <- read.csv("barplot_problem.csv")
library(ggplot2)
ggplot(dat, aes(x=counts.vc, y=duration, group=duration.type,
                colour=duration.type, fill=duration.type)) +
  stat_summary(fun.y=mean, geom="bar", position="stack") +
  stat_summary(fun.data=mean_cl_normal, geom="errorbar", colour="black") +
  facet_wrap(~speaker)


#https://www.jianshu.com/p/68efabd2e32b
#https://stackoverflow.com/questions/30872977/how-to-stack-error-bars-in-a-stacked-bar-plot-using-geom-errorbar



##2019.12.5

#https://mp.weixin.qq.com/s/l6B71b85vyYvIVsMn2eXjQ

##你需要堆叠柱状图添加bar吗？

##构造误差线的坐标ddply(df_res,"Attribute",transform,label_y = cumsum(Mean ))；
##第二个是重新排布堆叠柱状图不同分组因子水平，保证按照正确的方向填充：factor(df_res$Species,levels = c("virginica","versicolor","setosa"  ))。
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(plyr)


## 载入数据,这里是默认的鸢尾花数据
df <- iris
#数据宽边长
df <- melt(df, id="Species", variable.name="Attribute", value.name = "Size")
#设置出图颜色
mycol= brewer.pal(n = 12, name = "Set3")


## 数据统计均值、标准差、标准误
mean <- aggregate(df$Size, by=list(df$Species, df$Attribute), FUN=mean)
sd <- aggregate(df$Size, by=list(df$Species, df$Attribute), FUN=sd)
len <- aggregate(df$Size, by=list(df$Species, df$Attribute), FUN=length)
df_res <- data.frame(mean, sd=sd$x, len=len$x)
colnames(df_res) = c("Species", "Attribute", "Mean", "Sd", "Count")
head(df_res)
df_res$Se <- df_res$Sd/sqrt(df_res$Count) ### 计算标准差

#构造误差线坐标;?transform
df_res1 = ddply(df_res,"Attribute",transform,label_y = cumsum(Mean ));?ddply
head(df_res1)
#因子重新排列
df_res1$Species = factor(df_res$Species,levels = c("virginica","versicolor","setosa"  ))

### ggplot 绘图
ggplot(df_res1, aes(x=Attribute, y=Mean, fill=Species)) +
  geom_bar(stat="identity",color="black", width=.6) +
  scale_fill_manual(values = mycol) +
  geom_errorbar(aes(ymin=label_y-Sd, ymax=label_y +Sd), width=.2)
