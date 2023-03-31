setwd("C:/Users/24925/Desktop/Rdata")



#install.packages('survival')
##在导入临床数据之前在excel上把没有信息的都筛选掉，女0男1，都改成数字形式，列名改一下
comgene <- intersect(rownames(risk_all),rownames(phenotype_sur))#取交集
table(substr(comgene,14,16))
phenotype_sur <- cbind(phenotype_sur[comgene,],risk_all[comgene,]$riskScore)
colnames(phenotype_sur)[7] <- 'riskScore'
write.csv(phenotype_sur,file="phenotype_sur_score.csv")
rt <-phenotype_sur 


library(survival)
library(forestplot)#加载forestplot包
#样本名为第一列做行名，注意去除空白格，时间和生存状态分别放在第二第三列，后面为需要分析的变量
rt=read.table("duli.txt",header=T,sep="\t",check.names=F,row.names=1)
rt <- phenotype_sur_score1
#单因素cox
uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)

senlin <- read.table("uniCox.txt",header=T,sep="\t")#导入数据
senlin$`HRCI` <- ifelse(is.na(senlin$HR), "",
                           sprintf("%.2f (%.2f- %.2f)",
                                   senlin$HR, senlin$HR.95L, senlin$HR.95H))
senlin$pvalue <- sprintf("%.2f ",senlin$pvalue)


## 出现在森林图中的元素
tabletext <- cbind(c("Characteristics","\n",senlin$id),
                   c("HR(95%CI)","\n",senlin$HRCI),
                     c("P Value","\n",senlin$pvalue))

##绘制森林图
forestplot(labeltext=tabletext, 
           graph.pos=3,#为Pvalue箱线图所在的位置
           mean=c(NA,NA,senlin$HR),
           lower=c(NA,NA,senlin$HR.95L), upper=c(NA,NA,senlin$HR.95H),
           #定义标题
           title="Hazard Ratio Plot",
           ##定义x轴标题
           #xlab="    <---PCI Better---   ---Medical Therapy Better--->",
           #划线
           hrzl_lines=list("1" = gpar(lwd=1, col="black"),
             "2" = gpar(lwd=1, col="black"),
             "7" = gpar(lwd=1, col="black")),
           
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1.0),#字大小
                          ticks=gpar(cex=1.0),
                          xlab=gpar(cex = 1.0),
                          title=gpar(cex = 1.2)),
           ##fpColors函数设置颜色
           col=fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"),
           #箱线图中基准线的位置
           zero=1,
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           #箱子大小，线的宽度
           lwd.ci=1.8, boxsize=0.3,
           #箱线图两端添加小竖线，高度
           ci.vertices=TRUE, ci.vertices.height = 0.2)


#多因素cox
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)

senlin <- read.table("multiCox.txt",header=T,sep="\t")#导入数据
senlin$`HRCI` <- ifelse(is.na(senlin$HR), "",
                        sprintf("%.2f (%.2f- %.2f)",
                                senlin$HR, senlin$HR.95L, senlin$HR.95H))
senlin$pvalue <- sprintf("%.2f ",senlin$pvalue)


## 出现在森林图中的元素
tabletext <- cbind(c("Characteristics","\n",senlin$id),
                   c("HR(95%CI)","\n",senlin$HRCI),
                   c("P Value","\n",senlin$pvalue))

##绘制森林图
forestplot(labeltext=tabletext, 
           graph.pos=3,#为Pvalue箱线图所在的位置
           mean=c(NA,NA,senlin$HR),
           lower=c(NA,NA,senlin$HR.95L), upper=c(NA,NA,senlin$HR.95H),
           #定义标题
           title="Hazard Ratio Plot",
           ##定义x轴标题
           #xlab="    <---PCI Better---   ---Medical Therapy Better--->",
           #划线
           hrzl_lines=list("1" = gpar(lwd=1, col="black"),
                           "2" = gpar(lwd=1, col="black"),
                           "7" = gpar(lwd=1, col="black")),
           
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1.0),#字大小
                          ticks=gpar(cex=1.0),
                          xlab=gpar(cex = 1.0),
                          title=gpar(cex = 1.2)),
           ##fpColors函数设置颜色
           col=fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"),
           #箱线图中基准线的位置
           zero=1,
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           #箱子大小，线的宽度
           lwd.ci=1.8, boxsize=0.3,
           #箱线图两端添加小竖线，高度
           ci.vertices=TRUE, ci.vertices.height = 0.2)
