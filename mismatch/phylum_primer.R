###2019.5.14
#环形柱状图

#https://www.jianshu.com/p/410a909fe911
#https://www.r-graph-gallery.com/297-circular-barplot-with-groups/

#准备数据
df <- data.frame(individual=paste("Mister",seq(1,60),sep=""),value=sample(seq(10,100),60,replace=T))
df$id <- seq(1,nrow(df))
library(ggplot2)
#简易柱形图
p<-ggplot(df,aes(x=as.factor(id),y=value))+geom_bar(stat="identity",fill="blue")
#简易环状柱形图
p+coord_polar()

#环状图中间搞成空心，看起来好像美观一点
p+ylim(-100,120)+coord_polar()
#添加标签
p+coord_polar()+ylim(-100,120)+
  geom_text(aes(x=id,y=value+20,label=individual),size=3)+
  theme_minimal()+ylab("")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

#为参数hjust和angle赋予数据来调控标签的位置
df$angle <- 96-df$id*6
ggplot(df,aes(x=as.factor(id),y=value))+
  geom_bar(stat="identity",fill=alpha("blue",0.7))+
  coord_polar()+ylim(-100,120)+
  geom_text(aes(x=id,y=value+20,label=individual,angle=angle),
            size=3,hjust=0.2)+
  theme_minimal()+ylab("")+xlab("")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank())

#在完善一下
df$angle1<-ifelse(df$id<=30,96-df$id*6,96-df$id*6+180)
df$hjust<-ifelse(df$id<=30,0.2,1)
ggplot(df,aes(x=as.factor(id),y=value))+
  geom_bar(stat="identity",fill=alpha("blue",0.7))+
  coord_polar()+ylim(-100,120)+
  geom_text(aes(x=id,y=value+20,label=individual,
                angle=angle1,hjust=hjust),size=3)+
  theme_minimal()+ylab("")+xlab("")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank())

#value中添加缺失值，画图时就会出现空白的部分从而达到分割的目的
df1<-data.frame(individual=paste("Mister",seq(1,60),sep=""),
                value=rep(c(sample(60:100,9,replace=T),NA),6))
df1$id<-seq(1,nrow(df1))
df1
df1$angle<-df$angle1
df1$hjust<-df$hjust
df1
df1$fill<-c(rep("A",10),rep("B",10),rep("C",10),rep("D",10),rep("E",10),rep("F",10))
ggplot(df1,aes(x=as.factor(id),y=value))+
  geom_bar(stat="identity",aes(fill=fill))+
  coord_polar()+ylim(-100,120)+
  geom_text(aes(x=id,y=value+20,label=individual,
                angle=angle,hjust=hjust),size=3)+
  theme_minimal()+ylab("")+xlab("")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position="none")+
  scale_fill_manual(values=c("red","yellow","blue","green","orange","skyblue"))

#小知识点：ggplot2更改绘图区空白大小
#https://ggplot2.tidyverse.org/reference/element.html
theme(plot.margin=unit(c(1,1,1,1),'cm'))
#更改里面的数值即可
#比如可以比较一下以下两条命令的区别
df<-data.frame(A=1:10,B=10:1)
p<-ggplot(df,aes(x=A,y=B))+geom_point()
p+theme(plot.margin=unit(1,1,1,1),'cm')
p+theme(plot.margin=unit(2,2,2,2),'cm')






###############################
setwd("E:/桌面/ITS引物评价文章/2018.1.10 十二稿-Molcular ecology resources/文章修改文件/primer/")
x<-read.table(file="ITS1.txt",sep="\t",header=T)  
x
x$id<-seq(1,nrow(x))
n = 360/nrow(x)
m = round(nrow(x)/2)
x$angle<-ifelse(x$id<=m,96-x$id*n,96-x$id*n+180)  ##
x$hjust<-ifelse(x$id<=m,0.05,1) ##左右两部分的phylum名字距离
#x$fill<-c(rep("A",9),rep("B",2),rep("C",5)) #####每个引物分别调整三组的个数
x$fill<-c(rep("A",9),rep("B",2)) #####每个引物分别调整三组的个数
x
library(ggplot2)
p1 = ggplot(x,aes(x=as.factor(id),y=value))+
  geom_bar(stat="identity",aes(fill=fill))+
  coord_polar()+ylim(-80,180)+
  geom_text(aes(x=id,y=value+20,label=Phylum,
                angle=angle,hjust=hjust),size = 4,fontface = "bold",family = "serif")+  
  theme_minimal()+ylab("")+xlab("")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position="none")+
  scale_fill_manual(values=c("red","orange","skyblue"))
  #scale_fill_manual(values=c("red","orange"))
p1
# prepare a data frame for base lines
base_data=x %>% 
  group_by(fill) %>% 
  summarize(start=min(id)+1, end=max(id)) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
p1 = p1 + geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )
p1

#################
##优化，中间加数字和弧线。hjust一直有问题，不用了。
library(tidyverse)

# Create dataset
data<-read.table(file="ITS1ss.txt",sep="\t",header=T)  ;data

# Set a number of 'empty bar' to add at the end of each group
empty_bar=3
to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) = colnames(data)
to_add$group=rep(levels(data$group), each=empty_bar)
data=rbind(data, to_add)
data=data %>% arrange(group)
data$id=seq(1, nrow(data))
data
# Get the name and the y position of each label
label_data=data
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data=data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

# Make the plot
p = ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() ;p
p = p + geom_text(data=label_data, aes(x=id, y=value+10, label=Phylum, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

p

#############5.21
#累积柱状图
setwd("E:/桌面/ITS引物评价文章/2018.1.10 十二稿-Molcular ecology resources/primer/")
x<-read.table(file="ITS4-fun.txt",sep="\t",header=T)  
x
library(ggplot2)

#门默认按照字母顺序排列
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

p  = ggplot(x,aes(Phylum,percent,fill=type))+geom_bar(stat="identity",position="stack")
p = p + theme(axis.text= element_text(size=20, color="black", family  = "serif", face= "bold", vjust=0.5, hjust=0.5))+
        theme(axis.title = element_text(size=20, color="black", family  = "serif",face= "bold", vjust=0.5, hjust=0.5))
p = p + ylab("Percentage of covered phylum") + xlab("Phylum in database")+ 
        scale_y_continuous( breaks=seq(0,100,25))+
    guides(fill=guide_legend(title="ITS4-Fun")) + coord_flip()
p = p + plot_theme; p


####2019.8.30

##柱状图不够展示我的数据一气之下我把它掰弯了 

##https://mp.weixin.qq.com/s/XrNIpMUeV0rGbeZMW69hTA

##https://www.data-to-viz.com/graph/circularbarplot.html
##https://www.data-to-viz.com/graph/circularpacking.html