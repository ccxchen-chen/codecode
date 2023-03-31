
setwd("C:/Users/24925/Desktop/Rdata")

library(coin)
library(readxl)
library(openxlsx)
library(forestmodel)
library(survival)
library(forestplot)
rownames(lasso_gene) <- lasso_gene[,1]  

comgene <- intersect(rownames(lasso_gene),colnames(ucox_surv.expr_train))#取交集
rt<- cbind(ucox_surv.expr_train$OS,ucox_surv.expr_train$OS.time,ucox_surv.expr_train[,comgene])
colnames(rt)[1] <- 'fustat'#给生存状态和生存时间改个代码
colnames(rt)[2] <- 'futime'
data <- rt
#data <- read.xlsx("39lncrna.xlsx")
head(data)


rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)
rt[,"futime"]=rt[,"futime"]/365
cox <- coxph(Surv(futime, fustat) ~ ., data = rt)#test组数据整理需要这一步
cox=step(cox,direction = "both")#筛选基因
riskScore=risk_train$riskScore
summary=summary(cox)


#####计算mcox的HR值
colnames(summary$conf.int)
multi1<-as.data.frame(round(summary$conf.int[, c(1, 3, 4)], 1))
#一-3、multi2：提取：HR(95%CI)和P
library(tableone)
multi2<-ShowRegTable(cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#一-4.将两次提取结果合并成表；取名result
result <-cbind(multi1,multi2);result
#一-5.行名转为表格第一列，并给予命名"lncRNAs"
result<-tibble::rownames_to_column(result, var = "lncRNAs");result
write.csv(result,file="coxResultHR.csv")
#########


coxGene=rownames(summary$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))

write.csv(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
          file="risk_train.csv")
write.csv(cbind(id=coxGene,summary$coefficients),file="coxResult.csv")


####test组分成高低风险组####
rownames(coxResult) <- coxResult[,1]  
#rownames(ucox_surv.expr_test) <- ucox_surv.expr_test[,1]  
#ucox_surv.expr_test = ucox_surv.expr_test[,-1]
comgene <- intersect(rownames(coxResult),colnames(ucox_surv.expr_test))#取交集
rt_test<- cbind(ucox_surv.expr_test$OS,ucox_surv.expr_test$OS.time,ucox_surv.expr_test[,comgene])
colnames(rt_test)[1] <- 'fustat'#给生存状态和生存时间改个代码
colnames(rt_test)[2] <- 'futime'
rt_test[,"futime"]=rt_test[,"futime"]/365

#计算风险评分test
riskScore_test=risk_test$riskScore
risk=as.vector(ifelse(riskScore_test>median(riskScore_test),"high","low"))
risk_test <- cbind(id=rownames(cbind(rt_test[,outCol],riskScore_test,risk)),cbind(rt_test[,outCol],riskScore_test,risk))
colnames(risk_test)[12] <- 'riskScore'
write.csv(risk_test,file="risk_test.csv")
####all组分成高低风险组####

rownames(coxResult) <- coxResult[,1]  
#rownames(ucox_surv.expr_test) <- ucox_surv.expr_test[,1]  
#ucox_surv.expr_test = ucox_surv.expr_test[,-1]
comgene <- intersect(rownames(coxResult),colnames(ucox_surv.expr_01))#取交集
rt_all<- cbind(ucox_surv.expr_01$OS,ucox_surv.expr_01$OS.time,ucox_surv.expr_01[,comgene])
colnames(rt_all)[1] <- 'fustat'#给生存状态和生存时间改个代码
colnames(rt_all)[2] <- 'futime'
rt_all[,"futime"]=rt_all[,"futime"]/365
#计算风险评分all
riskScore_all=risk_all$riskScore
risk=as.vector(ifelse(riskScore_all>median(riskScore_all),"high","low"))
risk_all <- cbind(id=rownames(cbind(rt_all[,outCol],riskScore_all,risk)),cbind(rt_all[,outCol],riskScore_all,risk))
colnames(risk_all)[12] <- 'riskScore'
write.csv(risk_all,file="risk_all.csv")







#PCA数据分组#
exp_risk <- diff_gene_exp_fpkm
riskScore=predict(cox,type="risk",newdata=exp_risk)
risk=as.vector(ifelse(riskScore_test>median(riskScore),"high","low"))
table(rownames(exp_risk)==LINC01412)

