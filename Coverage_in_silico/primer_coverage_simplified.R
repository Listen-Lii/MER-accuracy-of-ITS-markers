
http://www.dxy.cn/bbs/thread/31834883#31834883
http://www.bio-info-trainee.com/926.html



# 输入引物输出长度，简并度及所有序列 -------------------------------------------------------

##直接读序列
dd = DNAStringSet("ANCTNMVYRWSKMDVHBN")  ##不识别U
rr = RNAStringSet("ANCUNMVYRWSKMDVHBN")  ##不识别T
##readDNAStringSet("") 打开一个文件
# 
# R: A or G
# Y: C or T
# S: G or C
# W: A or T
# K: G or T
# M: A or C
# B: C or G or T
# D: A or G or T
# H: A or C or T
# V: A or C or G
# N: any base

R= c("A,G")
Y= c("C","T")
S= c("G","C")
W= c("A","T")
K= c("G","T")
M= c("A","C")
B= c("C","G","T")
D= c("A","G","T")
H= c("A","C","T")
V= c("A","C","G")
N = c("A","G","C","T")

rm(list=ls())
library(Biostrings)
#primer = DNAStringSet("GTGARTCATCGARTCTTTG") 
primer = DNAString("GTGARTCATCGARTCTTTG") 
#nchar(primer)
length(primer)
for (i in 1:length(primer)) {
  i=4
  if ( as.character(primer[i])== ("A|T|C|G") ){
    i = 1+1
  }
  else{
    substr(primer[i]) ==
    
  }
}

as.character(primer[i])== (' "A"|"T" ')
A== ("A|T")


chartr("R","A",primer)
chartr("R","G",primer)
R= c("A,G")
Y= c("C","T")
S= c("G","C")
W= c("A","T")
K= c("G","T")
M= c("A","C")
B= c("C","G","T")
D= c("A","G","T")
H= c("A","C","T")
V= c("A","C","G")
N = c("A","G","C","T")
