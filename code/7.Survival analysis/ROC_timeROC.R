####ROC####

#准备R包
#install.packages("ROCR")
#install.packages("rms")
library(ROCR)
library(rms)
exp_sur <- risk_all

#构建ROC预测模型
ROC1 <- prediction(exp_sur$riskScore,exp_sur$fustat)   #构建ROC预测模型#SPP1可以改成别的基因，看这个基因对生存时间的预测能力 
ROC2 <- performance(ROC1,"tpr","fpr")   #计算预测模型的TPR/FPR值
AUC <- performance(ROC1,"auc")   #计算曲线下面积(AUC)值

AUC<- 0.775834 #改 根据结果对AUC进行赋值，从得出的结果里看

#1.4 绘制ROC曲线
plot(ROC2,
     col="red",   #曲线的颜色
     xlab="False positive rate", ylab="True positive rate",   #x轴和y轴的名称
     lty=1,lwd=3,
     main=paste("AUC=",AUC))
aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",ROC2$AUC),")"))
abline(0, 1, lty=2, lwd=3)   #绘制对角线
dev.off()

####timeROC####
#预测一年的生存率是多少，三年的生存率是多少
#setwd("timeROC")
#R包
#install.packages("timeROC")
#install.packages("survival")
library(timeROC)
library(survival)
library(tidyverse)


#2.3 构建ROC曲线函数
ROC3 <- timeROC(T=exp_sur$futime,   #结局时间
                delta=exp_sur$fustat,   #结局指标
                marker=exp_sur$riskScore,   #预测变量#使用的基因可以改变
                cause=1,   #阳性结局指标数值
                weighting="marginal",   #计算方法，默认为marginal
                times=c(1, 3, 5),   #时间点，选取1年，3年和5年的生存率
                iid=TRUE)
ROC3   #查看模型变量信息

#2.4 绘制ROC曲线
plot(ROC3,
     time=1, col="red",lwd=3,title=FALSE)   #time是时间点，col是线条颜色
plot(ROC3,
     time=3, col="green", add=TRUE,lwd=3,title=FALSE)   #add指是否添加在上一张图中
plot(ROC3,
     time=5, col="blue", add=TRUE,lwd=3,title=FALSE)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC3$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC3$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC3$AUC[3]))),
       col=c("green",'blue','red'),lwd=3,bty = 'n')
dev.off()
