library(survival)
library(rms)
library(ROCR)
library(Hmisc)
if (!require("survcomp")) {
  BiocManager::install("survcomp")
}
library(survcomp)
library(rms)
library(foreign)
library(survival)
com <- intersect(rownames(phenotype_sur_score),rownames(risk_all))#列名取交集，为了使两个表的列名变得相同
phenotype_sur_score$risk <- risk_all[com,]$risk

#读取数据
#TARGET<-read.table("TARGET_ver.txt",header=T,sep="\t")
TARGET <- phenotype_sur_score
#将数据转换成因子格式 
TARGET$age=as.numeric(TARGET$age)
#TARGET$age<-factor(TARGET$age,labels=c("<=65",">65"))
TARGET$Gender<-factor(TARGET$gender,labels=c("female","male"))
TARGET$stage<-factor(TARGET$stage,labels=c("stage 1","stage 2","stage 3","stage 4"))
#TARGET$T<-factor(TARGET$T,labels=c("T4","T3","T2","T1"))
#TARGET$N<-factor(TARGET$N,labels=c("N3","N2","N1","N0"))
#TARGET$M<-factor(TARGET$M,labels=c("M1","M0"))
TARGET$risk<-factor(TARGET$risk,labels=c("low","high"))

#1-year
cox1 <- cph(Surv(futime,fustat) ~ age  +gender+ stage ,surv=T,x=T, y=T,time.inc = 1*365*1,data=TARGET) 
#下面的m=120是设置分组的，看分成几个部分，注意不能被样本总数整除，要有余数
#B=1000，这个数认为越大越好，但是一般1000就够了，B过大，样本量过多，会让运算时间变长
cal <- calibrate(cox1, cmethod="KM", method="boot", u=1*365*1, m= 120, B=1000)

##c-index
coxph <- coxph(Surv(futime, fustat) ~ age + stage+gender  , data = TARGET)
coxph
cindex <- concordance.index(predict(coxph), surv.time = TARGET$futime, surv.event = TARGET$fustat)
cindex$c.index
## [1] 0.4335472
cindex$lower
## [1] 0.2439767
cindex$upper
## [1] 0.6447896


#画校准图
#pdf("calibrate有risk.pdf")
plot(cal,lwd=2,lty=1,errbar.col="black",xlim = c(0,1),ylim = c(0,1),xlab ="Nomogram-Predicted Probability of Survival",ylab="Actual Survival",col="green",sub=F)
#lines(cal,c('mean.predicted','KM'),type = 'a',lwd = 1,col ="black" ,pch = 18)
mtext(" ")
box(lwd = 0.5)
abline(0,1,lty = 3,lwd = 2,col = "black")




#3-year
cox1 <- cph(Surv(futime,fustat) ~ age + stage+gender ,surv=T,x=T, y=T,time.inc = 1*365*3,data=TARGET) 
#下面的m=120是设置分组的，看分成几个部分，注意不能被样本总数整除，要有余数
#B=1000，这个数认为越大越好，但是一般1000就够了，B过大，样本量过多，会让运算时间变长
cal <- calibrate(cox1, cmethod="KM", method="boot", u=1*365*3, m= 120, B=1000)


#画校准图
#pdf("calibrate3.pdf")
plot(cal,lwd=2,lty=1,errbar.col="black",xlim = c(0,1),ylim = c(0,1),xlab ="Nomogram-Predicted Probability of Survival",ylab="Actual Survival",col="blue",sub=F,add=TRUE)
#lines(cal,c('mean.predicted','KM'),type = 'a',lwd = 1,col ="black" ,pch = 18)
mtext(" ")
box(lwd = 0.5)
abline(0,1,lty = 3,lwd = 2,col = "black")


#5-year
cox1 <- cph(Surv(futime,fustat) ~ age + stage+gender ,surv=T,x=T, y=T,time.inc = 1*365*5,data=TARGET) 

cal <- calibrate(cox1, cmethod="KM", method="boot", u=1*365*5, m= 120, B=1000)


#画校准图
#pdf("calibrate5.pdf")
plot(cal,lwd=2,lty=1,errbar.col="black",xlim = c(0,1),ylim = c(0,1),xlab ="Nomogram-Predicted Probability of Survival",ylab="Actual Survival",col="red",sub=F,add=TRUE)
#lines(cal,c('mean.predicted','KM'),type = 'a',lwd = 1,col ="black" ,pch = 18)
mtext(" ")
box(lwd = 0.5)
abline(0,1,lty = 3,lwd = 2,col = "black")

legend("bottomright",
       c("Year-1", "Year-3", "Year-5"),
       col=c("green", "blue", "red"),
       lty=1, lwd=2)   #添加标签信息
legend("topleft",                                        #图例位置为上方
        legend="c-index=0.76502")

dev.off()
