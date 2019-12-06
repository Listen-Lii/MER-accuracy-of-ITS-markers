#读入文件
x<-read.table(file="mis.txt",sep="\t",header=T)  
x

#每组分别计算不同mismatch的标准偏差及坐标
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
#合并
da.ta = rbind(E_A,E_B,E_C,EK_A,EK_B,EK_C,S_A,S_B,S_C,SK_A,SK_B,SK_C)
data = cbind(facet=c(rep("Even community",18),rep("Staggered community",18)),da.ta)
data
data$level =rep(c("A","C","B"),12)

data2 = data %>% arrange(Type,level);data2

#画图
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
