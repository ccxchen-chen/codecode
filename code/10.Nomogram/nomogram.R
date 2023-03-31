####nomogram####
#导入数据，选取行名是样本名，列名是不同临床信息+生存时间、生存状态+风险分组、风险评分


comgene <- intersect(rownames(risk_all),rownames(phenotype_sur_score))#取交集
table(substr(comgene,14,16))
rt2 <- cbind(phenotype_sur_score[comgene,],risk_all[comgene,]$risk)
colnames(rt2)[6] <- 'risk'


#加载包
library(rms)
library(foreign)
library(survival)
###3. 设置参数
#意思是：现在rt2表格里的数据都是字符形式的，现在要把它变为因子型（相当于excel里的下拉框，要么是0，要么是1）
rt2$gender <- factor(rt2$gender,labels=c("F", "M"))
rt2$age=as.numeric(rt2$age)
#rt2$age <- factor(rt2$age,labels=c("<=65", ">65"))
rt2$stage <- factor(rt2$stage,labels=c("Stage1", "Stage2", "Stage3", "Stage4"))
#rt2$T <- factor(rt2$T,labels=c("T1", "T2", "T3", "T4"))
#rt2$M <- factor(rt2$M,labels=c("MX","M0", "M1"))
#rt2$N <- factor(rt2$N,labels=c("NX", "N0","N1"))#MX,NX不能留
rt2$risk <- factor(rt2$risk,labels=c("high","low"))
#rt2$race <- factor(rt2$race,labels=c("white", "black or african american","asian"))

#rt2 <- rt2[,-5]#这部是为了去掉没有预测价值的信息的列
rt2$futime <- rt2$futime/365


ddist <- datadist(rt2)
options(datadist='ddist')   #使用函数datadist()将数据打包

###4. 构建Cox回归模型
f <- cph(Surv(futime, fustat) ~age + stage +risk+gender, x=T, y=T, surv=T, data=rt2, time.inc=1)
surv <- Survival(f)
#需要去除没有预测价值的因素
###5. 构建Nomogram
nom2 <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(3, x), function(x) surv(5, x)), #数字1、2、3代表预测1、2、3年的生存数据
                 lp=F, funlabel=c("1-year survival", "3-year survival", "5-year survival"), 
                 maxscale=100, 
                 fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05))
plot(nom2)
dev.off()

