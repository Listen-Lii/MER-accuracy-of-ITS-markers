##2018.11
##计算引物在分类及序列上的覆盖度

###需要更新的地方
#1.R 本地blast
#2. 简并度和引物长度
#3. #,',和分列的问题

rm(list=ls(all=TRUE));gc()
setwd("E:/桌面/统计数据库物种")

# 匹配序列和物种信息--------------------------------------------------
# 从我们服务器上下载Qiime版本的silva数据库的silva_132_97.fasta这个没用。不清楚是涵盖哪些区域的序列。前两步不用做了。
# 这两步是刚开始做的时候最先用这个数据库写的代码。代码还是有用的。

#1.silva数据库,下载的物种文件,先将fungi单独挑出。后续用于统计总的分类信息

#y<-read.table(file="taxonomy_all_levels.txt",sep="\t",blank.lines.skip=F)  ##不加blank.lines.skip=F会报错，203行没有16个元素。blank.lines.skip 是否跳过空白行
#去掉所有的#和‘
y<-read.table(file="taxonomy_all_levels_without_().txt",sep="\t")          

#fungi = y[order(y[,5]),]  #第五列（fungi）排序
fungi <- subset(y, y[,5]=="D_3__Fungi") ;head(fungi,2) ##挑出所有的真菌 
#or fungi = y[y[,5]=="D_3__Fungi",]
#or fungi = y[which(y[,5]=="D_3__Fungi"), ]
#or library(dplyr); data <- filter(y, y[,5]=="D_3__Fungi")
write.table(fungi,file = "taxa_fungi.txt",sep="\t",row.names = F,col.names = F,quote=F) 
#quote=F输出内容不带引号。不加这个都会带引号。

###Silva_132_notes.txt
# raw_taxonomy.txt - these are the sequence IDs followed by the raw taxonomy strings directly pulled from the SILVA NR fasta file (will work with the -m blast assignment method, but not uclust/RDP)
# 
# taxonomy_7_levels.txt - This is the raw taxa, forced into exactly 7 levels as described in the preceding paragraph. This will work with all assignment methods
# 
# taxonomy_all_levels.txt - This is the raw taxa, expanded out to all levels present in any of the taxonomy strings (14 total levels). Will work with all assignment methods, but will use more memory than the 7 level taxonomy. Deeper levels of taxonomy, which will mostly come from Eukaryotes will require expansion of levels used with QIIME scripts, such as summarize_taxa.py.
# 
# consensus_taxonomy_7_levels.txt - This file is the same as the 7 levels, but uses the 100% consensus taxonomy (this is described in the “Consensus and Majority Taxonomies” section).
# 
# consensus_taxonomy_all_levels.txt - This file is the same as the all levels taxonomy, but uses the 100% consensus taxonomy (this is described in the “Consensus and Majority Taxonomies” section).
# 
# majority_taxonomy_7_levels.txt - This file is the same as the 7 levels, but uses the 90% majority taxonomy (this is described in the “Consensus and Majority Taxonomies” section).
# 
# majority_taxonomy_all_levels.txt - This file is the same as the all levels taxonomy, but uses the 90% majority taxonomy (this is described in the “Consensus and Majority Taxonomies” section).

#x<-read.table(file="majority_taxonomy_7_levels_tab.txt",sep="\t")
# D=unique(x[,5]);D   ##真菌那一级
# length(D)
# dim(D)=c(length(D),1)
# D = as.data.frame(D);D
# write.table(D,file = "taxa number.txt",sep="\t",col.names=NA )

#2.数据库fasta和物种信息匹配。物种信息作为fasta序列名称
library(Biostrings) 
s = readDNAStringSet("silva_132_97_18S.fna");head(s)
seqnames = s@ranges@NAMES

# taxa <-read.table(file="taxa_fungi.txt",sep="\t",row.names = 1)  
# taxaname = rownames(taxa)
taxaname = fungi[,1]

mat = match(taxaname,seqnames)
fungiseq = s[seqnames[mat]] ##从数据库中挑出有真菌物种信息的序列

##序列名称改为物种
#length(fungiseq)  ##3127
for (i in c(1:length(fungiseq))){
  fungiseq@ranges@NAMES[i]= paste(taxa[i,1],taxa[i,2],taxa[i,3],taxa[i,4],taxa[i,5],taxa[i,6],taxa[i,7],taxa[i,8],taxa[i,9],taxa[i,10],taxa[i,11],taxa[i,12],sep="\t")
}
head(fungiseq@ranges@NAMES,2)
max = max(fungiseq@ranges@width);max
writeXStringSet(fungiseq, file="fungiseq_silva.fasta",format="fasta",width=max+100) ##width 规定一行最大的字符数。注意输出后\r\n替换为\n
#?writeXStringSet


# Unite_coverage  ---------------------------------------------------------
#最新的数据库是2017.12
#unite数据库有四种。分别是带s(包含singletons)和不带s(不包含singletons)。
#developer是没有经过ITSx切除的，包含了部分LSU和SSU序列。我用的是developer带s的数据库。
#developer的存在对中间ITS引物区域的覆盖度几乎没有影响。而引物在LSU和SSU区间时又不需要用这个数据库来做。所以developer不是很重要。

#得到blastn结果后计算覆盖度
rm(list=ls(all=TRUE))
#setwd("E:/桌面/ITS引物覆盖度")
setwd("E:/桌面/统计数据库物种")
##总的序列数；在门及属水平上物种数。unite数据库对应物种信息
database1 = read.table(file="total_taxa_unite.txt",sep="\t")  ##6 is phylum, 10 is genus
seqnum = nrow(database);seqnum                    ##58048条序列  #length(unique(database[,3]))

##2019.5.24更新。仍有真菌之外的，先去掉。
database = database1[which(database1[,5]=="k__Fungi"),]
dphy=unique(database[,6])  ;length(dphy)          ####35个门  更新后20个门
dcla=unique(database[,7])  ;length(dcla)          #89
dord=unique(database[,8])  ;length(dord)          #249
dfam=unique(database[,9])  ;length(dfam)          #685
dgen=unique(database[,10]) ;length(dgen)          ###3626个属
dspe=unique(database[,11]) ;length(dspe)          #20543

phy.num = as.data.frame(table(database[,6]));phy.num
##注意，unidentified这种需要去掉。其实就是少了一个。所以在计算的时候加一个判断是否有unidentified就行了。
#数据库每一级必然有unidentified，计算时减一。
#balst结果要判断一下。这个其实对结果影响微乎其微。出现unidentified只能是低分类水平上。分子或分母减一对结果影响不大。
# unspe = which(grepl("^.*unidentified.*",database[,11],fixed=F))
# dspe=unique(database[,11][-unspe]) ;length(dspe)          #20542


##blast结果
setwd("E:/桌面/ITS引物覆盖度")
x<-read.table(file="ITS2.txt",sep="\t")  ##文件里的##注意替换掉##7 phylum, 11 genus
primerlength = 20;degeneracy = 0  ##长度20，简并0  2^0 =1 

x<-read.table(file="gITS7.txt",sep="\t")  ##文件里的;注意替换为tab ##7 phylum, 11 genus
primerlength = 19;degeneracy = 2  ##长度19，简并2  2^2=4

x<-read.table(file="5.8S.txt",sep="\t")  ##文件里的;注意替换为tab ##7 phylum, 11 genus
primerlength = 21;degeneracy = 5  ##长度21，简并5  2^5=32

# x<-read.table(file="E:/桌面/NRM-真菌综述/引物评价/gITS7ngs.txt",sep="\t")  ##文件里的;和|注意替换为tab ##7 phylum, 11 genus
# primerlength =19 ;degeneracy = 4  ##长度21，简并4  2^4=16

##mismatch = 0  筛选符合的序列  13列是identify，14列是aligh length,15列是mismatch
mismatch = 0
x_mis0 = subset(x,(x[,13]=="100")&(x[,14]==primerlength)&(x[,15]==mismatch))

phylum=unique(x_mis0[,7]) ;length(phylum) 
class=unique(x_mis0[,8]) ;length(class) 
order=unique(x_mis0[,9]) ;length(order) 
family=unique(x_mis0[,10]) ;length(family) 
genus=unique(x_mis0[,11]);length(genus)
species=unique(x_mis0[,12]) ;length(species) 
#table(unlist(x_mis0[,3]))  统计出现的次数
matchseq = length(unique(x_mis0[,3]));matchseq   ##注意，这个地方要取unique。否则有简并存在序列会重复计算。覆盖度可能会大于100%。
total_coverage = paste("total coverage is",matchseq*100/seqnum,"%");total_coverage  ##total coverage

#注意由于简并序列会多，因此先取出unique再统计匹配上的序列数
eachphy.num = as.data.frame(table(x_mis0[unique(x_mis0[,3]),][,7]));eachphy.num

##phylum
intphy = intersect(phylum,dphy)  ##交集
#phylum_coverage = paste("phylum coverage is",length(intphy)*100/length(dphy),"%");phylum_coverage  #门水平覆盖度
if ( length(grep("^.*unidentified.*",intphy,fixed=F))>0 ) { 
  phylum_coverage = paste("phylum coverage is",(length(intphy)-1)*100/(length(dphy)-1),"%") 
} else { phylum_coverage = paste("phylum coverage is",length(intphy)*100/(length(dphy)-1),"%") };phylum_coverage
intphy =c("intersectphylum number is", length(intphy),intphy)
difphy = setdiff(dphy,phylum)   ##取前一个里面不同的元素。注意顺序
difphy =c("diffphylum number is", length(difphy),difphy)
##class
intcla = intersect(class,dcla)  ##交集
#class_coverage = paste("class coverage is",length(intcla)*100/length(dcla),"%")
if ( length(grep("^.*unidentified.*",intcla,fixed=F))>0 ) { 
  class_coverage = paste("class coverage is",(length(intcla)-1)*100/(length(dcla)-1),"%") 
} else { class_coverage = paste("class coverage is",length(intcla)*100/(length(dcla)-1),"%") };class_coverage  
intcla =c("intersectclass number is", length(intcla),intcla)
difcla = setdiff(dcla,class)   ##取前一个里面不同的元素。注意顺序
difcla =c("diffclass number is", length(difcla),difcla)
##order
intord = intersect(order,dord)  ##交集
#order_coverage = paste("order coverage is",length(intord)*100/length(dord),"%");
if ( length(grep("^.*unidentified.*",intord,fixed=F))>0 ) { 
  order_coverage = paste("order coverage is",(length(intord)-1)*100/(length(dord)-1),"%") 
} else { order_coverage = paste("order coverage is",length(intord)*100/(length(dord)-1),"%") };order_coverage  
intord =c("intersectorder number is", length(intord),intord)
diford = setdiff(dord,order)   ##取前一个里面不同的元素。注意顺序
diford =c("difforder number is", length(diford),diford)
##family
intfam = intersect(family,dfam)  ##交集
#family_coverage = paste("family coverage is",length(intfam)*100/length(dfam),"%");
if ( length(grep("^.*unidentified.*",intfam,fixed=F))>0 ) { 
  family_coverage = paste("family coverage is",(length(intfam)-1)*100/(length(dfam)-1),"%") 
} else { family_coverage = paste("family coverage is",length(intfam)*100/(length(dfam)-1),"%") };family_coverage  
intfam =c("intersectfamily number is", length(intfam),intfam)
diffam = setdiff(dfam,family)   ##取前一个里面不同的元素。注意顺序
diffam =c("difffamily number is", length(diffam),diffam)
##genus
intgen = intersect(genus,dgen)  ##交集
#genus_coverage = paste("genus coverage is ",length(intgen)*100/length(dgen),"%")
if ( length(grep("^.*unidentified.*",intgen,fixed=F))>0 ) { 
  genus_coverage = paste("genus coverage is",(length(intgen)-1)*100/(length(dgen)-1),"%") 
} else { genus_coverage = paste("genus coverage is",length(intgen)*100/(length(dgen)-1),"%") };genus_coverage  #属水平覆盖度
intgen =c("intersectgenus number is", length(intgen),intgen)
difgen = setdiff(dgen,genus)   ##取前一个里面不同的元素。注意顺序
difgen =c("diffgenus number is", length(difgen),difgen)
##species
intspe = intersect(species,dspe)  ##交集
#species_coverage = paste("species coverage is",length(intspe)*100/length(dspe),"%");
if ( length(grep("^.*unidentified.*",intspe,fixed=F))>0 ) { 
  species_coverage = paste("species coverage is",(length(intspe)-1)*100/(length(dspe)-1),"%") 
} else { species_coverage = paste("species coverage is",length(intspe)*100/(length(dspe)-1),"%") };species_coverage 
intspe =c("intersectspecies number is", length(intspe),intspe)
difspe = setdiff(dspe,species)   ##取前一个里面不同的元素。注意顺序
difspe =c("diffspecies number is", length(difspe),difspe)

getwd()
options(max.print=100000000)
sink("5.8S_mismatch=0.txt")

# unite_mis0_out  ---------------------------------------------------------
total_coverage 
phylum_coverage 
class_coverage
order_coverage
family_coverage
genus_coverage
species_coverage
intphy
difphy
intcla
difcla
intord
diford
intfam
diffam
intgen
difgen
intspe
difspe
sink()

##mismatch = 1
mismatch = 1
x_mis1 = subset(x,((x[,14]==primerlength)&((x[,15]==mismatch)|(x[,15]==0)))
                
                |( (x[,14]==primerlength-1)&(x[,15]==0) ) )

#x_mis1 = subset(x,((x[,14]==primerlength)&((x[,15]== (mismatch|0) )))  这种写法就不行。会丢掉大部分。不清楚为什么
                
#                |( (x[,14]==primerlength-1)&(x[,15]==0) ) )    
#(x[,15]== (mismatch|0))  ==     ((x[,15]==mismatch)|(x[,15]==0))

#即
##align length = primer,mismatch=0  +
##aligh length = primer,mismtch=1   +
##align length = primer-1,mismtch = 0

phylum=unique(x_mis1[,7]) ;length(phylum) 
class=unique(x_mis1[,8]) ;length(class) 
order=unique(x_mis1[,9]) ;length(order) 
family=unique(x_mis1[,10]) ;length(family) 
genus=unique(x_mis1[,11]);length(genus)
species=unique(x_mis1[,12]) ;length(species) 
#matchseq = nrow(x_mis1);matchseq
matchseq = length(unique(x_mis1[,3]));matchseq   ##注意，这个地方要取unique。否则有简并存在序列会重复计算。覆盖度可能会大于100%。
total_coverage = paste("total coverage is",matchseq*100/seqnum,"%");total_coverage  ##total coverage

#注意由于简并序列会多，因此先取出unique再统计匹配上的序列数
eachphy.num = as.data.frame(table(x_mis1[unique(x_mis1[,3]),][,7]));eachphy.num

# mm = x_mis1[which(x_mis1[,7] == "p__GS01"),]
# length(unique(mm[,3]))

##phylum
intphy = intersect(phylum,dphy)  ##交集
#phylum_coverage = paste("phylum coverage is",length(intphy)*100/length(dphy),"%");phylum_coverage  #门水平覆盖度
if ( length(grep("^.*unidentified.*",intphy,fixed=F))>0 ) { 
  phylum_coverage = paste("phylum coverage is",(length(intphy)-1)*100/(length(dphy)-1),"%") 
} else { phylum_coverage = paste("phylum coverage is",length(intphy)*100/(length(dphy)-1),"%") };phylum_coverage
intphy =c("intersectphylum number is", length(intphy),intphy)
difphy = setdiff(dphy,phylum)   ##取前一个里面不同的元素。注意顺序
difphy =c("diffphylum number is", length(difphy),difphy)
##class
intcla = intersect(class,dcla)  ##交集
#class_coverage = paste("class coverage is",length(intcla)*100/length(dcla),"%")
if ( length(grep("^.*unidentified.*",intcla,fixed=F))>0 ) { 
  class_coverage = paste("class coverage is",(length(intcla)-1)*100/(length(dcla)-1),"%") 
} else { class_coverage = paste("class coverage is",length(intcla)*100/(length(dcla)-1),"%") };class_coverage  
intcla =c("intersectclass number is", length(intcla),intcla)
difcla = setdiff(dcla,class)   ##取前一个里面不同的元素。注意顺序
difcla =c("diffclass number is", length(difcla),difcla)
##order
intord = intersect(order,dord)  ##交集
#order_coverage = paste("order coverage is",length(intord)*100/length(dord),"%");
if ( length(grep("^.*unidentified.*",intord,fixed=F))>0 ) { 
  order_coverage = paste("order coverage is",(length(intord)-1)*100/(length(dord)-1),"%") 
} else { order_coverage = paste("order coverage is",length(intord)*100/(length(dord)-1),"%") };order_coverage  
intord =c("intersectorder number is", length(intord),intord)
diford = setdiff(dord,order)   ##取前一个里面不同的元素。注意顺序
diford =c("difforder number is", length(diford),diford)
##family
intfam = intersect(family,dfam)  ##交集
#family_coverage = paste("family coverage is",length(intfam)*100/length(dfam),"%");
if ( length(grep("^.*unidentified.*",intfam,fixed=F))>0 ) { 
  family_coverage = paste("family coverage is",(length(intfam)-1)*100/(length(dfam)-1),"%") 
} else { family_coverage = paste("family coverage is",length(intfam)*100/(length(dfam)-1),"%") };family_coverage  
intfam =c("intersectfamily number is", length(intfam),intfam)
diffam = setdiff(dfam,family)   ##取前一个里面不同的元素。注意顺序
diffam =c("difffamily number is", length(diffam),diffam)
##genus
intgen = intersect(genus,dgen)  ##交集
#genus_coverage = paste("genus coverage is ",length(intgen)*100/length(dgen),"%")
if ( length(grep("^.*unidentified.*",intgen,fixed=F))>0 ) { 
  genus_coverage = paste("genus coverage is",(length(intgen)-1)*100/(length(dgen)-1),"%") 
} else { genus_coverage = paste("genus coverage is",length(intgen)*100/(length(dgen)-1),"%") };genus_coverage  #属水平覆盖度
intgen =c("intersectgenus number is", length(intgen),intgen)
difgen = setdiff(dgen,genus)   ##取前一个里面不同的元素。注意顺序
difgen =c("diffgenus number is", length(difgen),difgen)
##species
intspe = intersect(species,dspe)  ##交集
#species_coverage = paste("species coverage is",length(intspe)*100/length(dspe),"%");
if ( length(grep("^.*unidentified.*",intspe,fixed=F))>0 ) { 
  species_coverage = paste("species coverage is",(length(intspe)-1)*100/(length(dspe)-1),"%") 
} else { species_coverage = paste("species coverage is",length(intspe)*100/(length(dspe)-1),"%") };species_coverage 
intspe =c("intersectspecies number is", length(intspe),intspe)
difspe = setdiff(dspe,species)   ##取前一个里面不同的元素。注意顺序
difspe =c("diffspecies number is", length(difspe),difspe)

getwd()
options(max.print=100000000)
sink("5.8S_mismatch=1.txt")

# unite_mis1_out ----------------------------------------------------------
total_coverage 
phylum_coverage 
class_coverage
order_coverage
family_coverage
genus_coverage
species_coverage
intphy
difphy
intcla
difcla
intord
diford
intfam
diffam
intgen
difgen
intspe
difspe
sink()



# LSU_coverage ------------------------------------------------------------
rm(list=ls())
setwd("E:/桌面/统计数据库物种")

###LSU和SSU用的都是Silva 132版本数据库。最新版本是2017.12.13
###一定注意这里是RNA序列。包含有U。如果读DNA会把所有U全部丢弃。
###Silva_132_notes.txt中提到convert U characters to T characters

library(Biostrings) 
##LSU = readDNAStringSet("Galaxy42-[SILVA_132_LSURef_tax_silva_trunc.fasta].fasta")
LSU = readRNAStringSet("Galaxy42-[SILVA_132_LSURef_tax_silva_trunc.fasta].fasta");head(LSU) 
seqnames = LSU@ranges@NAMES ; head(seqnames)
length(seqnames)  ##198843

#zz = gsub(";","\t",seqnames);head(zz)  ##;替换为\t
##fixed=F 表示为正则表达式
length(grep("^.*Fungi.*",zz,fixed=F))  ##5129   ##从数据库中挑出有真菌物种信息的序列数
#LSU_fungi = LSU[seqnames[grep("^.*Fungi.*",zz,fixed=F)]]  ##fungiseq = s[seqnames[mat]]

###从数据库中挑出有真菌物种信息的序列数
LSU_fungi = LSU[seqnames[grep("^.*Fungi.*",seqnames,fixed=F)]]  ##fungiseq = s[seqnames[mat]]
 # LSU_fungi@ranges@NAMES = gsub(";","\t",LSU_fungi@ranges@NAMES)
 # LSU_fungi@ranges@NAMES = gsub(" ","\t",LSU_fungi@ranges@NAMES)  ##;和空格都替换为\t
 
##output taxa information
LSU_taxa = gsub(";","\t",LSU_fungi@ranges@NAMES) 
print(paste("含有#和‘的个数为",length(grep("'|#",LSU_fungi@ranges@NAMES))))
LSU_taxa = gsub("#|'","",LSU_taxa)
LSU_taxa = as.data.frame(LSU_taxa)
write.table(LSU_taxa,"LSU_taxa.txt",sep="\t",col.names = F,row.names=F,quote=F)
#write.table(LSU_taxa,"LSU_taxa_Primerprostor.txt",sep="\t",col.names = F,row.names=F,quote=F)

#excel中将第一列按照空格分列。
##将物种信息对齐


##空格匹配为，
LSU_fungi@ranges@NAMES = gsub(" ",",",LSU_fungi@ranges@NAMES)  
head(LSU_fungi@ranges@NAMES);length(LSU_fungi)  ##5129  没问题
max = max(LSU_fungi@ranges@width)  ##最长3159
print(paste("含有#和‘的个数为",length(grep("'|#",LSU_fungi@ranges@NAMES))))
LSU_fungi@ranges@NAMES = gsub("#|'","",LSU_fungi@ranges@NAMES)  ##替换掉#和’
writeXStringSet(LSU_fungi, file="LSU_fungi.fasta",format="fasta",width=max+1000) ##width 规定一行最大的字符数。
##注意输出后\r\n替换为\n
getwd()
##输出文件做blast即可。

##blast结果先替换，为空格，；为tab.  ##特殊字符都去掉#()
##excel第一列分列并手动补齐列再导入文件
#结果也需要对齐
##LSU for ITS4R and ITS4-FunR
rm(list=ls())
setwd("E:/桌面/ITS引物覆盖度")

#x<-read.table(file="LSU_ITS4R.txt",sep="\t",fill=T,col.names=paste("V",1:22,sep=""))  ##3456
x<-read.table(file="LSU_ITS4R_all_unidentified.txt",sep="\t")  ##3456
primerlength = 20;degeneracy = 0  

# x<-read.table(file="LSU_ITS4Fun_R.txt",sep="\t",fill=T,col.names=paste("V",1:22,sep="")) ##6972
# primerlength = 27;degeneracy = 1  

###LSU_ITS4Fun去掉前两个碱基
x<-read.table(file="Galaxy101-[blastn.txt]_ITS4Fun_without first 2 base_all_unidentified.txt",sep="\t",fill=T,col.names=paste("V",1:22,sep="")) ##6972
primerlength = 25;degeneracy = 1  

##总的序列数；在门水平上物种数.LSU目后面就是种，没有科和属。
#LSU_taxa = read.table(file = "LSU_taxa.txt",sep = "\t")  ##直接读显示不完全
#LSU for ITS4R and ITS4Fun-R
setwd("E:/桌面/统计数据库物种")
LSU_taxa = read.table(file = "LSU_taxa_all_unidentified.txt",sep = "\t") ##必须去掉#‘这两个特殊符号才能显示完全
seqnum = nrow(LSU_taxa);seqnum                    ##5129条序列  
dphy=unique(LSU_taxa[,7])  ;length(dphy)          ####9个门
dcla=unique(LSU_taxa[,9])  ;length(dcla)          #25
dord=unique(LSU_taxa[,10])  ;length(dord)         #60
dspe=unique(LSU_taxa[,11])  ;length(dspe)         #1889

phy.num = as.data.frame(table(LSU_taxa[,7]));phy.num

##mismatch = 0  筛选符合的序列  13列是identify，14列是aligh length,15列是mismatch
mismatch = 0
x_mis0 = subset(x,(x[,13]=="100")&(x[,14]==primerlength)&(x[,15]==mismatch))

phylum=unique(x_mis0[,8]) ;length(phylum)   ##8 is phylum
class=unique(x_mis0[,10]) ;length(class)
order=unique(x_mis0[,11]) ;length(order)
species=unique(x_mis0[,12]) ;length(species)
#genus=unique(x_mis0[,11]);length(genus)
#matchseq = nrow(x_mis0);matchseq
###table(unlist(ss[,3]))  统计出现的次数
matchseq = length(unique(x_mis0[,2]));matchseq   ##注意，这个地方要取unique。否则有简并存在序列会重复计算。覆盖度可能会大于100%。
total_coverage = paste("total coverage is",matchseq*100/seqnum,"%");total_coverage  ##total coverage

####将unidentified”, “uncultured”, “metagenome” “eukaryotic picoplankton environmental sample” “fungal sp.”, and “incertae sedis”
####这些全部替换为cultured即可。这样少了很多判断，计算更简单。

#注意由于简并序列会多，因此先取出unique再统计匹配上的序列数
eachphy.num = as.data.frame(table(x_mis0[unique(x_mis0[,2]),][,8]));eachphy.num

##phylum
intphy = intersect(phylum,dphy)  ##交集
phylum_coverage = paste("phylum coverage is",(length(intphy)-1)*100/(length(dphy)-1),"%");phylum_coverage  #门水平覆盖度
intphy =c("intersectphylum number is", length(intphy),intphy)
difphy = setdiff(dphy,phylum)   ##取前一个里面不同的元素。注意顺序
difphy =c("diffphylum number is", length(difphy),difphy)
##class
intcla = intersect(class,dcla)  ##交集
class_coverage = paste("class coverage is",(length(intcla)-1)*100/(length(dcla)-1),"%");class_coverage  
intcla =c("intersectclass number is", length(intcla),intcla)
difcla = setdiff(dcla,class)   ##取前一个里面不同的元素。注意顺序
difcla =c("diffclass number is", length(difcla),difcla)
##order
intord = intersect(order,dord)  ##交集
order_coverage = paste("order coverage is",(length(intord)-1)*100/(length(dord)-1),"%");order_coverage  
intord =c("intersectorder number is", length(intord),intord)
diford = setdiff(dord,order)   ##取前一个里面不同的元素。注意顺序
diford =c("difforder number is", length(diford),diford)
##species
intspe = intersect(species,dspe)  ##交集
species_coverage = paste("species coverage is",(length(intspe)-1)*100/(length(dspe)-1),"%");species_coverage 
intspe =c("intersectspecies number is", length(intspe),intspe)
difspe = setdiff(dspe,species)   ##取前一个里面不同的元素。注意顺序
difspe =c("diffspecies number is", length(difspe),difspe)
##genus
# intgen = intersect(genus,dgen)  ##交集
# genus_coverage = paste("genus coverage is ",(length(intgen)-1)*100/(length(dgen)-1),"%");genus_coverage  #属水平覆盖度
# intgen =c("intersectgenus number is", length(intgen),intgen)
# difgen = setdiff(dgen,genus)   ##取前一个里面不同的元素。注意顺序
# difgen =c("diffgenus number is", length(difgen),difgen)

getwd()
options(max.print=100000000)
sink("ITS4_mismatch=0.txt")
total_coverage 
phylum_coverage 
class_coverage
order_coverage
species_coverage
intphy
difphy
intcla
difcla
intord
diford
intspe
difspe
#intgen
#difgen
sink()

##mismatch = 1
mismatch = 1
x_mis1 = subset(x,((x[,14]==primerlength)&((x[,15]==mismatch)|(x[,15]==0)))
                
                |( (x[,14]==primerlength-1)&(x[,15]==0) ) )

#x_mis1 = subset(x,((x[,14]==primerlength)&((x[,15]== (mismatch|0) )))  这种写法就不行。会丢掉大部分。不清楚为什么

#                |( (x[,14]==primerlength-1)&(x[,15]==0) ) )    
#(x[,15]== (mismatch|0))  ==     ((x[,15]==mismatch)|(x[,15]==0))

#即
##align length = primer,mismatch=0  +
##aligh length = primer,mismtch=1   +
##align length = primer-1,mismtch = 0

phylum=unique(x_mis1[,8]) ;length(phylum) 
class=unique(x_mis1[,10]) ;length(class)
order=unique(x_mis1[,11]) ;length(order)
species=unique(x_mis1[,12]) ;length(species)
#genus=unique(x_mis1[,11]);length(genus)
#matchseq = nrow(x_mis1);matchseq
matchseq = length(unique(x_mis1[,2]));matchseq   ##注意，这个地方要取unique。否则有简并存在序列会重复计算。覆盖度可能会大于100%。
total_coverage = paste("total coverage is",matchseq*100/seqnum,"%");total_coverage  ##total coverage

#注意由于简并序列会多，因此先取出unique再统计匹配上的序列数
eachphy.num = as.data.frame(table(x_mis1[unique(x_mis1[,2]),][,8]));eachphy.num

##phylum
intphy = intersect(phylum,dphy)  ##交集
phylum_coverage = paste("phylum coverage is",(length(intphy)-1)*100/(length(dphy)-1),"%");phylum_coverage  #门水平覆盖度
intphy =c("intersectphylum number is", length(intphy),intphy)
difphy = setdiff(dphy,phylum)   ##取前一个里面不同的元素。注意顺序
difphy =c("diffphylum number is", length(difphy),difphy)
##class
intcla = intersect(class,dcla)  ##交集
class_coverage = paste("class coverage is",(length(intcla)-1)*100/(length(dcla)-1),"%");class_coverage  
intcla =c("intersectclass number is", length(intcla),intcla)
difcla = setdiff(dcla,class)   ##取前一个里面不同的元素。注意顺序
difcla =c("diffclass number is", length(difcla),difcla)
##order
intord = intersect(order,dord)  ##交集
order_coverage = paste("order coverage is",(length(intord)-1)*100/(length(dord)-1),"%");order_coverage  
intord =c("intersectorder number is", length(intord),intord)
diford = setdiff(dord,order)   ##取前一个里面不同的元素。注意顺序
diford =c("difforder number is", length(diford),diford)
##species
intspe = intersect(species,dspe)  ##交集
species_coverage = paste("species coverage is",(length(intspe)-1)*100/(length(dspe)-1),"%");species_coverage 
intspe =c("intersectspecies number is", length(intspe),intspe)
difspe = setdiff(dspe,species)   ##取前一个里面不同的元素。注意顺序
difspe =c("diffspecies number is", length(difspe),difspe)
##genus
# intgen = intersect(genus,dgen)  ##交集
# intgen =c("intersectgenus number is", length(intgen),intgen)
# genus_coverage = paste("genus coverage is ",length(intgen)*100/length(dgen),"%");genus_coverage  #属水平覆盖度
# difgen = setdiff(dgen,genus)   ##取前一个里面不同的元素。注意顺序
# difgen =c("diffgenus number is", length(difgen),difgen)

getwd()
options(max.print=100000000)
sink("ITS4_mismatch=1.txt")
total_coverage 
phylum_coverage 
class_coverage
order_coverage
species_coverage
#genus_coverage
intphy
difphy
intcla
difcla
intord
diford
intspe
difspe
#intgen
#difgen
sink()



# SSU_coverage ----------------------------------------------------------
#SSU类似
rm(list=ls())
setwd("E:/桌面/统计数据库物种")

library(Biostrings) 
#SSU = readDNAStringSet("SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta")
SSU = readRNAStringSet("SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta");head(SSU)
seqnames = SSU@ranges@NAMES ; head(seqnames)
length(seqnames)  ##695171

#zz = gsub(";","\t",seqnames);head(zz)  ##;替换为\t
length(grep("^.*Fungi.*",seqnames,fixed=F))  ##15748   ##从数据库中挑出有真菌物种信息的序列
#length(grep("^.*Fungi.*",zz,fixed=F))  ##15748   ##从数据库中挑出有真菌物种信息的序列
#SSU_fungi = SSU[seqnames[grep("^.*Fungi.*",zz,fixed=F)]]  ##fungiseq = s[seqnames[mat]]
SSU_fungi = SSU[seqnames[grep("^.*Fungi.*",seqnames,fixed=F)]]  ##fungiseq = s[seqnames[mat]]

##output taxa information
SSU_taxa = gsub(";","\t",SSU_fungi@ranges@NAMES) 
print(paste("含有#和‘的个数为",length(grep("'|#",SSU_fungi@ranges@NAMES))))
SSU_taxa = gsub("#|'","",SSU_taxa)
SSU_taxa = as.data.frame(SSU_taxa)
write.table(SSU_taxa,"SSU_taxa_new.txt",sep="\t",col.names = F,row.names=F,quote=F)
#excel中将第一列按照空格分列。
##将物种信息对齐

# SSU_fungi@ranges@NAMES = gsub(";","\t",SSU_fungi@ranges@NAMES)
# SSU_fungi@ranges@NAMES = gsub(" ","\t",SSU_fungi@ranges@NAMES)  ##;和空格都替换为\t
SSU_fungi@ranges@NAMES = gsub(" ",",",SSU_fungi@ranges@NAMES)   ##空格匹配为，
head(SSU_fungi@ranges@NAMES);length(SSU_fungi)  ##15748
print(paste("含有#和‘的个数为",length(grep("'|#",SSU_fungi@ranges@NAMES))))
SSU_fungi@ranges@NAMES = gsub("#|'","",SSU_fungi@ranges@NAMES)  ##替换掉#和’
max = max(SSU_fungi@ranges@width)  ##最长2227
writeXStringSet(SSU_fungi, file="SSU_fungi.fasta",format="fasta",width=max+1000) 
##注意输出后\r\n替换为\n
##输出文件运行blast

##blast结果先替换，为空格，；为tab.  ##特殊字符都去掉#()
##excel第一列分列并手动补齐列再导入文件

#LSU for ITS1F
rm(list=ls())
setwd("E:/桌面/ITS引物覆盖度")
#setwd("E:/桌面/统计数据库物种")
x<-read.table(file="ITS1F_all_unidentified.txt",sep="\t")  
primerlength = 22;degeneracy = 0  


SSU_taxa = read.table(file = "SSU_taxa_new_all_unidentified.txt",sep = "\t") ##必须去掉#‘这两个特殊符号才能显示完全
seqnum = nrow(SSU_taxa);seqnum                    ##15746条序列  
dphy=unique(SSU_taxa[,6])  ;length(dphy)        ####17个门
dcla=unique(SSU_taxa[,8])  ;length(dcla)        ##61
dord=unique(SSU_taxa[,9])  ;length(dord)        ##179
dspe=unique(SSU_taxa[,12]) ;length(dspe)        ##6179个种

phy.num = as.data.frame(table(SSU_taxa[,6]));phy.num

##mismatch = 0  筛选符合的序列  15列是identify，16列是aligh length,17列是mismatch
mismatch = 0
x_mis0 = subset(x,(x[,15]=="100")&(x[,16]==primerlength)&(x[,17]==mismatch))

phylum=unique(x_mis0[,8]) ;length(phylum)   ##8 is phylum  12
class=unique(x_mis0[,10]) ;length(class)   #47
order=unique(x_mis0[,11]) ;length(order)    #117
species=unique(x_mis0[,14]) ;length(species)  #2363
#genus=unique(x_mis0[,13]);length(genus)
matchseq = nrow(x_mis0);matchseq
###table(unlist(ss[,3]))  统计出现的次数
matchseq = length(unique(x_mis0[,2]));matchseq   ##注意，这个地方要取unique。否则有简并存在序列会重复计算。覆盖度可能会大于100%。
total_coverage = paste("total coverage is",matchseq*100/seqnum,"%");total_coverage  ##total coverage

#注意由于简并序列会多，因此先取出unique再统计匹配上的序列数
eachphy.num = as.data.frame(table(x_mis0[unique(x_mis0[,2]),][,8]));eachphy.num

##phylum
intphy = intersect(phylum,dphy)  ##交集
phylum_coverage = paste("phylum coverage is",(length(intphy)-1)*100/(length(dphy)-1),"%");phylum_coverage  #门水平覆盖度
intphy =c("intersectphylum number is", length(intphy),intphy)
difphy = setdiff(dphy,phylum)   ##取前一个里面不同的元素。注意顺序
difphy =c("diffphylum number is", length(difphy),difphy)
##class
intcla = intersect(class,dcla)  ##交集
class_coverage = paste("class coverage is",(length(intcla)-1)*100/(length(dcla)-1),"%");class_coverage  
intcla =c("intersectclass number is", length(intcla),intcla)
difcla = setdiff(dcla,class)   ##取前一个里面不同的元素。注意顺序
difcla =c("diffclass number is", length(difcla),difcla)
##order
intord = intersect(order,dord)  ##交集
order_coverage = paste("order coverage is",(length(intord)-1)*100/(length(dord)-1),"%");order_coverage  
intord =c("intersectorder number is", length(intord),intord)
diford = setdiff(dord,order)   ##取前一个里面不同的元素。注意顺序
diford =c("difforder number is", length(diford),diford)
##species
intspe = intersect(species,dspe)  ##交集
species_coverage = paste("species coverage is",(length(intspe)-1)*100/(length(dspe)-1),"%");species_coverage 
intspe =c("intersectspecies number is", length(intspe),intspe)
difspe = setdiff(dspe,species)   ##取前一个里面不同的元素。注意顺序
difspe =c("diffspecies number is", length(difspe),difspe)
##genus
# intgen = intersect(genus,dgen)  ##交集
# intgen =c("intersectgenus number is", length(intgen),intgen)
# genus_coverage = paste("genus coverage is ",length(intgen)*100/length(dgen),"%");genus_coverage  #属水平覆盖度
# difgen = setdiff(dgen,genus)   ##取前一个里面不同的元素。注意顺序
# difgen =c("diffgenus number is", length(difgen),difgen)

getwd()
options(max.print=100000000)
sink("SSU_ITSOF-T=0.txt")
total_coverage 
phylum_coverage 
class_coverage
order_coverage
species_coverage
#genus_coverage
intphy
difphy
intcla
difcla
intord
diford
intspe
difspe
#intgen
#difgen
sink()



##mismatch = 1
mismatch = 1
x_mis1 = subset(x,((x[,16]==primerlength)&((x[,17]==mismatch)|(x[,17]==0)))
                
                |( (x[,16]==primerlength-1)&(x[,17]==0) ) )

#x_mis1 = subset(x,((x[,14]==primerlength)&((x[,15]== (mismatch|0) )))  这种写法就不行。会丢掉大部分。不清楚为什么

#                |( (x[,14]==primerlength-1)&(x[,15]==0) ) )    
#(x[,15]== (mismatch|0))  ==     ((x[,15]==mismatch)|(x[,15]==0))

#即
##align length = primer,mismatch=0  +
##aligh length = primer,mismtch=1   +
##align length = primer-1,mismtch = 0

phylum=unique(x_mis1[,8]) ;length(phylum)   ##8 is phylum  14
class=unique(x_mis1[,10]) ;length(class)   #54
order=unique(x_mis1[,11]) ;length(order)    #152
species=unique(x_mis1[,14]) ;length(species)  #3630
matchseq = nrow(x_mis1);matchseq
###table(unlist(ss[,3]))  统计出现的次数
matchseq = length(unique(x_mis1[,2]));matchseq   ##注意，这个地方要取unique。否则有简并存在序列会重复计算。覆盖度可能会大于100%。
total_coverage = paste("total coverage is",matchseq*100/seqnum,"%");total_coverage  ##total coverage

#注意由于简并序列会多，因此先取出unique再统计匹配上的序列数
eachphy.num = as.data.frame(table(x_mis1[unique(x_mis1[,2]),][,8]));eachphy.num

##phylum
intphy = intersect(phylum,dphy)  ##交集
phylum_coverage = paste("phylum coverage is",(length(intphy)-1)*100/(length(dphy)-1),"%");phylum_coverage  #门水平覆盖度
intphy =c("intersectphylum number is", length(intphy),intphy)
difphy = setdiff(dphy,phylum)   ##取前一个里面不同的元素。注意顺序
difphy =c("diffphylum number is", length(difphy),difphy)
##class
intcla = intersect(class,dcla)  ##交集
class_coverage = paste("class coverage is",(length(intcla)-1)*100/(length(dcla)-1),"%");class_coverage  
intcla =c("intersectclass number is", length(intcla),intcla)
difcla = setdiff(dcla,class)   ##取前一个里面不同的元素。注意顺序
difcla =c("diffclass number is", length(difcla),difcla)
##order
intord = intersect(order,dord)  ##交集
order_coverage = paste("order coverage is",(length(intord)-1)*100/(length(dord)-1),"%");order_coverage  
intord =c("intersectorder number is", length(intord),intord)
diford = setdiff(dord,order)   ##取前一个里面不同的元素。注意顺序
diford =c("difforder number is", length(diford),diford)
##species
intspe = intersect(species,dspe)  ##交集
species_coverage = paste("species coverage is",(length(intspe)-1)*100/(length(dspe)-1),"%");species_coverage 
intspe =c("intersectspecies number is", length(intspe),intspe)
difspe = setdiff(dspe,species)   ##取前一个里面不同的元素。注意顺序
difspe =c("diffspecies number is", length(difspe),difspe)
##genus
# intgen = intersect(genus,dgen)  ##交集
# intgen =c("intersectgenus number is", length(intgen),intgen)
# genus_coverage = paste("genus coverage is ",length(intgen)*100/length(dgen),"%");genus_coverage  #属水平覆盖度
# difgen = setdiff(dgen,genus)   ##取前一个里面不同的元素。注意顺序
# difgen =c("diffgenus number is", length(difgen),difgen)

getwd()
options(max.print=100000000)
sink("SSU_ITS9MUNngs_mismatch=1.txt")
total_coverage 
phylum_coverage 
class_coverage
order_coverage
species_coverage
#genus_coverage
intphy
difphy
intcla
difcla
intord
diford
intspe
difspe
#intgen
#difgen
sink()




# 序列替换 --------------------------------------------------------------------
#Biostrings::dna2rna(LSU)    ##rna2dna,cDNA这些函数已经没有了
#cDNA(LSU)

# chartr("U", "T", LSU)    ##U替换为T   # ?chartr
# qqq = sub("U","T",LSU);head(qqq)
#这两种替换方式都会改变数据结构

# 分割数据 --------------------------------------------------------------------
# library(stringr)
# seqqq = str_split_fixed(names,"\t",n = Inf)
# ?str_split_fixed

# 读取长度不同的数据 ---------------------------------------------------------------
##其实这个用处不大。读进来数据也会错位。
##我是在excel中排序后生动改的。

#####################读取不同长度的数据
rm(list=ls())
setwd("E:/桌面/ITS引物覆盖度/")
test = read.table(file = "test_file.txt",sep="\t")  ##直接读会报错
test = read.table(file = "test_file.txt",sep="\t",fill=T) ##这样会错行，只根据前5行判断
test
##如果文件少于5行，根据全文判断长度。超过5行只根据前5行判断长度。
#http://bbs.pinggu.org/thread-2137520-1-1.html
#这种方法只能针对比较少的列，人为找到最长的进行规定
test = read.table(file = "test_file.txt",sep="\t",fill=T,col.names=paste("col",1:9,sep="")) 
test
##如果有几百万列，不知道哪一列最长呢？
#先不加fill参数，从报错中找到line 1 did not have 18 elements。以这个值开始进行循环。
#依次加一，如果不行，则会有一些行在某个值之后全是空值和NA。
#空值和NA都有即为限制条件

test = read.table(file = "test_file.txt",sep="\t") 
##输入的文件没有空值和NA
##line 和 element根据报错来
line = 1
element = 5
for (i in line+1:nrow(test)){
  i=5
  if(length(grep("^\\S$",test[i,],fixed=F)) & length(grep("^\\S*NA$",test[i,],fixed=F))>0 )
  {
    element =element + 1
    test = read.table(file = "test_file.txt",sep="\t",fill=T,col.names=paste("V",1:element,sep="")) 
  }
  i = i+1
}

element
test = read.table(file = "test_file.txt",sep="\t",fill=T,col.names=paste("V",1:element,sep="")) 

# grepl("^\\S$",its4r[24,],fixed=F)      ##是否有空值  ##grep返回下标，grepl返回逻辑值
# 
# grepl("^\\S*NA$",its4r[24,],fixed=F)   ##是否有NA


grepl("^\\s$",its4r[5,],fixed=F)    #问题是有东西的行也会判断有空值

grep("\\s",its4r[24,3],fixed=F)
class(its4r)
as.vector(its4r)



# R lical blast -----------------------------------------------------------
##R 本地blast
rm(list=ls())
setwd("E:/培训ppt/blast/bin")
#setwd("E:/桌面/统计数据库物种")
getwd()

system("E:/培训ppt/blast/bin/blastn.exe -query E:/桌面/ITS引物覆盖度/other_primer.txt -db E:/桌面/统计数据库物种/fungiseq_silva.fasta -out outblast.txt")

#建立数据库
##特殊符号要去掉，如×
#在我们网站上是看不到这些warming的。所以会少一些序列。
system("makeblastdb.exe -in E:/桌面/统计数据库物种/sh_unite.fasta -dbtype nucl -out E:/桌面/统计数据库物种/Unite")
#blastn
##文件要改成fasta的格式。引物那一行加>
#只能写在一行里。多行会报错
outblast = system("blastn.exe -query E:/桌面/ITS引物覆盖度/unite_primer.txt -db E:/桌面/统计数据库物种/Unite -out E:/桌面/统计数据库物种/outblast.txt -outfmt 7 -task blastn -max_target_seqs 100000000")
outblast

##blast不会自动判断简并
##必须提前手动判断
