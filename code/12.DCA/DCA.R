install.packages('ggDCA')

com <- intersect(rownames(phenotype_sur_score),rownames(risk_all))#列名取交集，为了使两个表的列名变得相同
phenotype_sur_score$risk <- risk_all[com,]$risk

com <- intersect(rownames(phenotype_sur_score),rownames(TCGA_KIRC_phenotype_sur))#列名取交集，为了使两个表的列名变得相同
phenotype_sur_score$age <- TCGA_KIRC_phenotype_sur[com,]$age
write.csv(phenotype_sur_score,file = "phenotype_sur_score2.csv")





#setwd("C:/Users/24925/Desktop/Rdata")

library(ggDCA)
library(rms)
library(foreign)

#risk_score需要设123等级
#训练集
data1<- phenotype_sur_score4
#data1<-read.table("allnom.txt",header=T,sep="\t")
bc <- na.omit(data1)#删除缺失值

#验证集
data2<-read.table("controlclai.txt",header=T,sep="\t")
bc <- na.omit(data2)

#####生成3个模型
ddist <- datadist(data1)
options(datadist='ddist')
age<-cph(Surv(futime,fustat)~age,bc)
#riskScore<-cph(Surv(futime,fustat)~riskScore,bc)
stage<-cph(Surv(futime,fustat)~stage,bc)
risk<-cph(Surv(futime,fustat)~risk,bc)
gender<-cph(Surv(futime,fustat)~gender,bc)
#race<-cph(Surv(futime,fustat)~race,bc)
#T<-cph(Surv(futime,fustat)~T,bc)
#M<-cph(Surv(futime,fustat)~M,bc)
#N<-cph(Surv(futime,fustat)~N,bc)

d_train <- dca(age,stage,risk,gender,times=60)####多个模型5年后生存率
ggplot(d_train)


