#提取3Dgene
setwd("D:/work/MF2021gepeng/pancancer_rawdata/3Dgenerawdata2")
library(oligo)
library(preprocessCore)
rm(list = ls())
path <- "D:/work/MF2021gepeng/pancancer_rawdata/3Dgenerawdata2/GSE139031_RAW580"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ paste(path,x,sep='/')})  #sapply应用完函数后返回矩阵
data <- lapply(filePath, function(x){read.table(x,sep="\t",skip = 7,header = TRUE)[,c(5,7,9)]}) #lapply应用完函数后返回列表
e <- matrix(c(1:(length(fileNames)*length(data[[1]][,1]))),ncol =length(fileNames))
for(i in 1:length(fileNames)){e[,i] <- data[[i]][,2]-data[[i]][,3]}
row.names(e) <- data[[1]][,1]
# q<- e[grep("Negative Control*",row.names(e)),]#阴性对照探针
# for(i in  1:length(fileNames)){a <- sort(q[,i]);l=a[0.05*length(a)];u=a[0.95*length(a)];result <- q[,i][q[,i]>l&q[,i]<=u];dd <- mean(result)+2*sd(result);e[,i][e[,i]<=dd] <-0;e[,i] <- e[,i]-mean(result) }

for(i in 1:length(fileNames)){e[,i][e[,i]<=0] <- min(e[,i][e[,i]>0])-0.1}
#c <- log2(e)
#which(is.na(c))
#d=normalize.quantiles(c)
#med.d2 <- medpolish(d)
#h <- d-med.d2[["residuals"]]
#row.names(h) <- row.names(e)
#summary(medp.d2[,1])
#for(i in 1:lengh(fileNames)){for (j in 0:31){for(k in 1:100){if(e[100*j+k,i]<q[j+1,i]){e[100*j+k,i] <-q[j+1,i]}}}}#将每列小于QC值换成QC值
#for(i in 1:length(fileNames)){for (j in 0:15){for(k in 1:200){if(e[200*j+k,i]<q[j+1,i]){e[200*j+k,i] <-q[j+1,i]}}}}#将每列小于QC值换成QC值
#which(is.na(e))#查NA
h <- basicRMA(e,row.names(e),TRUE,FALSE)
which(is.na(h))#查NA



g <- h[grep("hsa",row.names(h)),]#只保留hsa的
g <- round(g,3)
summary(g[,1])
which(is.na(g))
for(i in 1:length(fileNames)){g[,i][is.na(g[,i])] <- min(na.omit(g[,i]))-0.05}
m <- read.csv('D:/work/MF2021gepeng/pancancer_rawdata/3Dgenerawdata15/qcGSE139031.2565.580.csv',header = F,row.names = 1)
group <- m[1,]
colnames(g) <- group
#一次就行
# qcb3dgene <- read.csv('D:/work/MF2021gepeng/pancancer_rawdata/3Dgenerawdata15/qcb3Dgene.csv')
# rowname <- qcb3dgene[,1]
# rowname <- as.character(rowname)

g1<- g[match(rowname,row.names(g)),]
dim(g1)

colnames(g1) <- group
write.csv(g1,"GSE139031_RAW580.csv")
write.csv(group,"GSE139031_RAW580group.csv")




n <- m[1,1:dim(m)[2]]
table(as.character(n))

name <- GSE73002_eSet@phenoData@data[["Group"]]
table(name)
reduce <- which(name=="Prostate Disease")
h <- h[,-reduce]
dim(h)

colnames(h) <- n
write.csv(h,"qcGSE139031.2565.580.csv")
#for(i in 1:length(fileNames)){g[,i][is.na(g[,i])] <- 2}#除NA
#med.h <- medpolish(h)
#medp.h <- h-med.h[["residuals"]]
#which(is.na(g))找na值
#sum(is.na(g))计算na总数

#library(data.table)
library(miRBaseVersions.db)
head(select(miRBaseVersions.db,keys = row.names(M)[5],keytype = "MIMAT",columns = c("ACCESSION","NAME","VERSION")),1)
#medpolish.example
deaths <-
  +     rbind(c(14,15,14),
              +           c( 7, 4, 7),
              +           c( 8, 2,10),
              +           c(15, 9,10),
              +           c( 0, 2, 0))

dimnames(deaths) <- list(c("1-24", "25-74", "75-199", "200++", "NA"),
                         +                          paste(1973:1975))
med.d <- medpolish(deaths)
d2 <- h[,which(n=="Benign Ovarian Diseases")]
med.d2 <- medpolish(d2)
medp.d2 <- d2-med.d2[["residuals"]]
h <- read.csv("GSE59856.csv",row.names = 1)

m <- read.csv("D:/work/MF2021gepeng/pancancer_rawdata/36 expression matrices and disease conditions/3Dgenerawdata15GSE59856.csv",header=FALSE,row.names = 1)

d1 <- read.csv("GSE59856_RAW571.csv",row.names = 1,header = FALSE)
d2 <- read.csv("GSE73002_RAW4113.csv",row.names = 1,header = FALSE)
d3 <- read.csv("GSE85677_RAW125.csv",row.names = 1,header = FALSE)
d4 <- read.csv("GSE85679_RAW281.csv",row.names = 1,header = FALSE)
d5 <- read.csv("GSE106817_RAW4046.csv",row.names = 1,header = FALSE)
d6 <- read.csv("GSE110651_RAW147.csv",row.names = 1,header = FALSE)
d7 <- read.csv("GSE112264_RAW1591.csv",row.names = 1,header = FALSE)
d8 <- read.csv("GSE113486_RAW972.csv",row.names = 1,header = FALSE)
d9 <- read.csv("GSE119159_RAW139.csv",row.names = 1,header = FALSE)
d10 <- read.csv("GSE119892_RAW66.csv",row.names = 1,header = FALSE)
d11 <- read.csv("GSE122497_RAW5531.csv",row.names = 1,header = FALSE)
d12 <- read.csv("GSE124158_RAW1412.csv",row.names = 1,header = FALSE)
d13 <- read.csv("GSE134108_RAW79.csv",row.names = 1,header = FALSE)
d14 <- read.csv("GSE137140_RAW3924.csv",row.names = 1,header = FALSE)
d15 <- read.csv("GSE139031_RAW580.csv",row.names = 1,header = FALSE)
g1 <- d1[1,]
g2 <- d2[1,]
g3 <- d3[1,]
g4 <- d4[1,]
g5 <- d5[1,]
g6 <- d6[1,]
g7 <- d7[1,]
g8 <- d8[1,]
g9 <- d9[1,]
g10 <- d10[1,]
g11 <- d11[1,]
g12 <- d12[1,]
g13<- d13[1,]
g14<- d14[1,]
g15<- d15[1,]
d1 <- d1[-1,]
d2 <- d2[-1,]
d3 <- d3[-1,]
d4 <- d4[-1,]
d5 <- d5[-1,]
d6 <- d6[-1,]
d7 <- d7[-1,]
d8 <- d8[-1,]
d9 <- d9[-1,]
d10 <- d10[-1,]
d11 <- d11[-1,]
d12 <- d12[-1,]
d13 <- d13[-1,]
d14 <- d14[-1,]
d15 <- d15[-1,]
rowname <- row.names(d10)[na.omit(match(row.names(d1),row.names(d10)))]
my_merge <- function(df1, df2){ cbind(df1, df2)}
c1 <- d1[match(rowname,row.names(d1)),]
c2 <- d2[match(rowname,row.names(d2)),]
c3 <- d3[match(rowname,row.names(d3)),]
c4 <- d4[match(rowname,row.names(d4)),]
c5 <- d5[match(rowname,row.names(d5)),]
c6 <- d6[match(rowname,row.names(d6)),]
c7 <- d7[match(rowname,row.names(d7)),]
c8 <- d8[match(rowname,row.names(d8)),]
c9 <- d9[match(rowname,row.names(d9)),]
c10<- d10[match(rowname,row.names(d10)),]
c11 <- d11[match(rowname,row.names(d11)),]
c12 <- d12[match(rowname,row.names(d12)),]
c13 <- d13[match(rowname,row.names(d13)),]
c14 <- d14[match(rowname,row.names(d14)),]
c15 <- d15[match(rowname,row.names(d15)),]
data <- list(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15)
b<- Reduce(my_merge, data)
group <- c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15)
group <- as.character(group)
colnames(b) <- group
v <- colnames(b)
table(v)
for(i in 1:length(v)){if(v[i]==" Benign Bone and Soft Tissue Tumor"){v[i]<- "Benign Bone and Soft Tissue Tumor"}}
for(i in 1:length(v)){if(v[i]=="Benign Breast Disease"){v[i]<- "Benign Breast Diseases"}}
for(i in 1:length(v)){if(v[i]=="Benign Prostatic Diseases"){v[i]<- "Benign Prostate Diseases"}}
for(i in 1:length(v)){if(v[i]=="No  Intracrainial Malignant Tumors (Infarction)"){v[i]<- "Intracrainial Benign Diseases"}}
for(i in 1:length(v)){if(v[i]==" Intracrainial Malignant Tumors"){v[i]<- "Intracrainial Malignant Tumors"}}
table(v)
colnames(b) <- v
write.csv(b,"newqcb3Dgene.csv")
t <- v
for(i in 1:length(t)){if(t[i]=="Benign Bone and Soft Tissue Tumor"){t[i]<- "Non-cancer Diseases"}}
for(i in 1:length(t)){if(t[i]=="Benign Breast Diseases"){t[i]<- "Non-cancer Diseases"}}
for(i in 1:length(t)){if(t[i]=="Benign Ovarian Diseases"){t[i]<- "Non-cancer Diseases"}}
for(i in 1:length(t)){if(t[i]=="Benign Pancreatic or Biliary Tract Diseases"){t[i]<- "Non-cancer Diseases"}}
for(i in 1:length(t)){if(t[i]=="Benign Prostate Diseases"){t[i]<- "Non-cancer Diseases"}}
for(i in 1:length(t)){if(t[i]=="Chronic Hepatitis C"){t[i]<- "Non-cancer Diseases"}}
for(i in 1:length(t)){if(t[i]=="Intracrainial Benign Diseases"){t[i]<- "Non-cancer Diseases"}}
for(i in 1:length(t)){if(t[i]!="Healthy"&t[i]!="Non-cancer Diseases"){t[i]<- "Tumors"}}
table(t)
t <- t(t)
write.csv(t,"qcb3Dgene3kinds.csv")
q <- t
for(i in 1:length(q)){if(q[i]=="Non-cancer Diseases"){q[i]<- "Healthy"}}
table(q)
write.csv(q,"qcb3Dgene2kinds.csv")
s <- v




tb <- t(b)
write.csv(tb,"tb3Dgene.csv")
bn=normalize.quantiles(as.matrix(b))
for(i in 1:23484){bn[,i] <- round(as.numeric(bn[,i]),3)}
write.csv(bn,"qcbn3Dgene.csv")
for(i in 1:23484){b[,i] <- round(as.numeric(b[,i]),3)}


#新测试集与训练集做分位数归一化，顺序混乱
setwd('D:/work/MF2021gepeng/pancancer_rawdata/验证3Dgene数据集')
d1 <- read.csv("g1.csv",row.names = 1,header = FALSE)
d2 <- read.csv("qcb3Dgene.csv",row.names = 1,header = FALSE)
g1 <- d1[1,]
g2 <- d2[1,]
d1 <- d1[-1,]
d2 <- d2[-1,]
rowname <- row.names(d2)[na.omit(match(row.names(d1),row.names(d2)))]
my_merge <- function(df1, df2){ cbind(df1, df2)}
c1 <- d1[match(rowname,row.names(d1)),]
c2 <- d2[match(rowname,row.names(d2)),]
data <- list(c1,c2)
b<- Reduce(my_merge, data)
group <- c(g1,g2)
group <- as.character(group)
colnames(b) <- group
v <- colnames(b)
table(v)
library(oligo)
datab=as.data.frame(lapply(b,as.numeric))
h <- basicRMA(as.matrix(datab),row.names(datab),TRUE,FALSE)
row.names(h) <- row.names(b)
colnames(h) <- colnames(b)


which(is.na(b))#查NA

h=normalize.quantiles(as.matrix(datab))
ceshi <- h[,23485:26418]
colnames(ceshi) <- g1
row.names(ceshi) <- row.names(c1)
write.csv(ceshi,"ceshi.csv")



# which(is.na(h))#查NA
# for(i in 1:135){h[,i] <- round(as.numeric(h[,i]),3)}
# med.h <- medpolish(h)
# medp.h <- h-med.h[["residuals"]]
# write.csv(medp.h,"dn.csv")

#med.d2 <- medpolish(d)
#h <- d-med.d2[["residuals"]]


setwd("D:/work/MF2021gepeng/pancancer_rawdata/验证3Dgene数据集")
library(oligo)
library(preprocessCore)
rm(list = ls())
path <- "D:/work/MF2021gepeng/pancancer_rawdata/验证3Dgene数据集/GSE113740"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ paste(path,x,sep='/')})  #sapply应用完函数后返回矩阵
data <- lapply(filePath, function(x){read.table(x,sep="\t",skip = 7,header = TRUE)[,c(5,7,9)]}) #lapply应用完函数后返回列表
e <- matrix(c(1:(length(fileNames)*length(data[[1]][,1]))),ncol =length(fileNames))
for(i in 1:length(fileNames)){e[,i] <- data[[i]][,2]-data[[i]][,3]}
row.names(e) <- data[[1]][,1]
q<- e[grep("Negative Control*",row.names(e)),]#阴性对照探针
for(i in  1:length(fileNames)){a <- sort(q[,i]);l=a[0.05*length(a)];u=a[0.95*length(a)];result <- q[,i][q[,i]>l&q[,i]<=u];dd <- mean(result)+2*sd(result);e[,i][e[,i]<=dd] <-0;e[,i] <- e[,i]-mean(result) }
for(i in 1:length(fileNames)){e[,i][e[,i]<=0] <- min(e[,i][e[,i]>0])-0.1}
h <- basicRMA(e,row.names(e),TRUE,FALSE)
which(is.na(h))#查NA
g <- h[grep("hsa",row.names(h)),]#只保留hsa的
g <- round(g,3)
summary(g[,1])
m <- read.table("GSE113740.txt",sep='\t',header = FALSE)
group <- m
for(i in 1:length(group)){if(grepl("Gastric Cancer",group[i])){group[i]<- "Gastric Cancer"}}
for(i in 1:length(group)){if(grepl("Non-cancer",group[i])){group[i]<- "Non-cancer"}}
for(i in 1:length(group)){if(grepl("Colorectal Cancer",group[i])){group[i]<- "Colorectal Cancer"}}
for(i in 1:length(group)){if(grepl("Esophageal Cancer",group[i])){group[i]<- "Esophageal Cancer"}}
group <- as.character(group)
table(as.character(group))
qcb3dgene <- read.csv('D:/work/MF2021gepeng/pancancer_rawdata/3Dgenerawdata15/qcb3Dgene.csv')
rowname <- qcb3dgene[,1]
rowname <- as.character(rowname)
g1<- g[match(rowname,row.names(g)),]
dim(g1)
colnames(g1) <- group
write.csv(g1,"g1GSE113740.csv")
write.csv(group,"group.csv")


# setwd("D:/work/MF2021gepeng/pancancer_rawdata/153Dgene")
# a1 <-read.csv("qcb3Dgene.csv",row.names = 1,header = FALSE) 
# a2 <-read.csv("g1GSE113740.csv",row.names = 1,header = FALSE) 
# g1 <- a1[1,]
# g2 <- a2[1,]
# d1 <- a1[-1,]
# d2 <- a2[-1,]
# rowname <- row.names(d1)
# d <- cbind(d1,d2)
# group <- c(g1,g2)
# group <- as.character(group)
# colnames(d) <- group
# v <- colnames(d)
# table(v)
# library(preprocessCore)
# datab=as.data.frame(lapply(d,as.numeric))
# h=normalize.quantiles(as.matrix(datab))
# row.names(h) <- row.names(d)
# colnames(h) <- colnames(d)
# g3GSE113740 <- h[,23485:25301]
# write.csv(g3GSE113740,"g3GSE113740.csv")
# qcqc3Dgene <- h[,1:23484]
# write.csv(qcqc3Dgene,"qcqc3Dgene.csv")
# 
# #qcqcb
# rm(list=ls())
# library(preprocessCore)
# a1 <-read.csv("qcb3Dgene.csv",row.names = 1,header = FALSE) 
# a2 <-read.csv("g1GSE113740.csv",row.names = 1,header = FALSE) 
# g1 <- a1[1,]
# g2 <- a2[1,]
# d1 <- a1[-1,]
# d2 <- a2[-1,]
# rowname <- row.names(d1)
# datab1=as.data.frame(lapply(d1,as.numeric))
# h=normalize.quantiles(as.matrix(datab1))
# row.names(h) <- row.names(d1)
# write.csv(h,"qcqcb3Dgene.csv")
# d <- cbind(h,d2)
# group <- c(g1,g2)
# group <- as.character(group)
# datab=as.data.frame(lapply(d,as.numeric))
# h=normalize.quantiles(as.matrix(datab))
# row.names(h) <- row.names(d)
# colnames(h) <- group
# g3GSE113740 <- h[,23485:25301]
# write.csv(g3GSE113740,"qcqcb113740.csv")

# #qc测试集
# rm(list=ls())
# library(preprocessCore)
# a1 <-read.csv("qcqcb3Dgene.csv",row.names = 1,header = FALSE) 
# a2 <-read.csv("g1.csv",row.names = 1,header = FALSE) 
# g1 <- a1[1,]
# g2 <- a2[1,]
# d1 <- a1[-1,]
# d2 <- a2[-1,]
# rowname <- row.names(d1)
# d <- cbind(d1,d2)
# group <- c(g1,g2)
# group <- as.character(group)
# colnames(d) <- group
# datab=as.data.frame(lapply(d,as.numeric))
# h=normalize.quantiles(as.matrix(datab))
# row.names(h) <- row.names(d)
# colnames(h) <- colnames(d)
# g3GSE164174 <- h[,23485:dim(h)[2]]
# write.csv(g3GSE164174,"qcqcb164174.csv")
# 
# #qc其他芯片
# rm(list=ls())
# setwd("D:/work/MF2021gepeng/pancancer_rawdata/153Dgene")
# library(preprocessCore)
# a1 <-read.csv("qcb3Dgene.csv",row.names = 1,header = FALSE) 
# a2 <-read.csv("bagilent3.csv",row.names = 1,header = FALSE) 
# g1 <- a1[1,]
# g2 <- a2[1,]
# d1 <- a1[-1,]
# d2 <- a2[-1,]
# rowname <- row.names(d2)
# for (i in 1:length(rowname)){
#   if (rowname[i]=="hsa-miR-199a-3p"){rowname[i] <-"hsa-miR-199a-3p, hsa-miR-199b-3p"  }
#   if (rowname[i]=="hsa-miR-365a-3p"){rowname[i] <-"hsa-miR-365a-3p, hsa-miR-365b-3p"  }
#   if (rowname[i]=="hsa-miR-3689a-5p"){rowname[i] <-"hsa-miR-3689a-5p, hsa-miR-3689b-5p, hsa-miR-3689e"  }
#   if (rowname[i]=="hsa-miR-3689b-3p"){rowname[i] <-"hsa-miR-3689b-3p, hsa-miR-3689c"  }
#   if (rowname[i]=="hsa-miR-516a-3p"){rowname[i] <-"hsa-miR-516b-3p, hsa-miR-516a-3p"   }
#   if (rowname[i]=="hsa-miR-517a-3p"){rowname[i] <-"hsa-miR-517a-3p, hsa-miR-517b-3p"  }
#   if (rowname[i]=="hsa-miR-518a-5p"){rowname[i] <-"hsa-miR-527, hsa-miR-518a-5p"  }
#   if (rowname[i]=="hsa-miR-518e-5p"){rowname[i] <-"hsa-miR-519c-5p, hsa-miR-523-5p, hsa-miR-518e-5p, hsa-miR-522-5p, hsa-miR-519a-5p, hsa-miR-519b-5p"  }
#   if (rowname[i]=="hsa-miR-548aa"){rowname[i] <-"hsa-miR-548aa, hsa-miR-548t-3p"   }
#   if (rowname[i]=="hsa-miR-548ai"){rowname[i] <-"hsa-miR-548ai, hsa-miR-570-5p"  }
#   if (rowname[i]=="hsa-miR-548aj-5p"){rowname[i] <-"hsa-miR-548g-5p, hsa-miR-548x-5p, hsa-miR-548aj-5p"  }
#   if (rowname[i]=="hsa-miR-548am-5p"){rowname[i] <-"hsa-miR-548c-5p, hsa-miR-548o-5p, hsa-miR-548am-5p"  }
#   if (rowname[i]=="hsa-miR-548h-3p"){rowname[i] <-"hsa-miR-548z, hsa-miR-548h-3p"   }
# }
# row.names(d2) <- rowname
# rowname0 <- row.names(d2)[na.omit(match(row.names(d1),row.names(d2)))]
# c1 <- d1[match(rowname0,row.names(d1)),]
# c1=as.data.frame(lapply(c1,as.numeric))
# c2 <- d2[match(rowname0,row.names(d2)),]
# c2=as.data.frame(lapply(c2,as.numeric))
# row.names(c1) <- rowname
# row.names(c2) <- rowname
# group1 <- as.character(g1)
# group2 <- as.character(g2)
# colnames(c1) <- group1
# colnames(c2) <- group2
# write.csv(c1,"rawqcb3Dgene.csv")
# write.csv(c2,"rawqcbagilent.csv")

# d <- cbind(c1,c2)
# group <- c(g1,g2)
# group <- as.character(group)
# datab=as.data.frame(lapply(d,as.numeric))
# gc()
# h=normalize.quantiles(as.matrix(datab))
# row.names(h) <- row.names(d)
# colnames(h) <- group
# h <- round(h,2)
# write.csv(h,"qcbagilent3agilent3d3.csv")
# algroups <-colnames(h) 
# agilent3d <- h[,1:23484]
# qcbagilent <- h[,23485:dim(h)[2]]
# write.csv(agilent3d,"agilent3d3.csv")
# write.csv(qcbagilent,"qcbagilent3.csv")
# which(is.na(agilent3d))
# which(is.na(c1))
# which(is.na(qcbagilent))
# 
# rm(list=ls())
# 
# setwd("D:/work/MF2021gepeng/pancancer_rawdata/153Dgene")
# a <- read.csv("qcbagilent3agilent3d3.csv",header = FALSE,row.names = 1)

# rm(list=ls())
# a <- read.csv("rawqcb3Dgene.csv",row.names = 1,header = FALSE)
# b <- read.csv("rawqcbagilent.csv",row.names = 1,header = FALSE)
# group1 <- a[1,]
# group1 <- as.character(group1)
# group2 <- b[1,]
# group2 <- as.character(group2)
# a=a[-1,]
# a=as.data.frame(lapply(a,as.numeric))
# b=b[-1,]
# b=as.data.frame(lapply(b,as.numeric))
# ahealthy <- a[,c(199:205)]
# ahealthy <- ahealthy[,c(1,2)]
# colnames(ahealthy) <- rep("Healthy",2)
# colnames(bhealthy) <- rep("Healthy",2)
# bhealthy <- b[,c(1:7)]
# bhealthy <- bhealthy[,c(1,2)]
# boxplot(c(ahealthy,bhealthy))
# which(group1=="Healthy")
# which(group2=="Lung Cancer")
# bluncancer <- b[,c(1558,1560)]
# colnames(aluncancer) <- rep("Lung Cancer",2)
# colnames(bluncancer) <- rep("Lung Cancer",2)
# colnames(bluncancer) <- group2[c(1558,1560)]
# which(group1=="Lung Cancer")
# aluncancer <- a[,c(5512,5514)]
# boxplot(c(ahealthy,bhealthy,aluncancer,bluncancer))
# normalize1 <- function(x){return((x-min(x))/(max(x)-min(x)))}
# # ncahealthy <- as.data.frame(lapply(ahealthy, normalize1))
# 
# 
# 
# rm(list=ls())
# a <- read.csv("newqcb3Dgene3kinds.csv",header = T,row.names = 1)
# a <- as.character(a)
# a <- t(a)
# table(a)
# for(i in 1:length(a)){if(a[i]=="Healthy"){a[i] <- 0}}
# for(i in 1:length(a)){if(a[i]=="Non-cancer Diseases"){a[i] <- 1}}
# for(i in 1:length(a)){if(a[i]=="Tumors"){a[i] <- 2}}
# write.csv(a,"newqcb3Dgene3kindsn.csv",)
# 
# rm(list=ls())
# a <- read.csv("newqcb3Dgene.csv",header = F,row.names = 1)
# a <- a[1,]
# a <- as.character(a)
# a <- t(a)
# table(a)
# for(i in 1:length(a)){if(a[i]=="Healthy"){a[i] <- 0}}
# for(i in 1:length(a)){if(a[i]=="Non-cancer Diseases"){a[i] <- 1}}
# for(i in 1:length(a)){if(a[i]=="Tumors"){a[i] <- 2}}
# write.csv(a,"newqcb3Dgene3kindsn.csv",)






