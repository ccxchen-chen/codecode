####LASSO COX回归模型的构建####
setwd("lasso")

#install.packages("glmnet")
#install.packages("survival")
library("glmnet")
library("survival")

topgene <- cox_topgene
ucox_surv.expr <- ucox_surv.expr_train
rownames(topgene) <- topgene[,1]  
topgene = topgene[,-1]
comgene <- intersect(rownames(topgene),colnames(ucox_surv.expr))#取交集
rt<- cbind(ucox_surv.expr$OS,ucox_surv.expr$OS.time,ucox_surv.expr[,comgene])
colnames(rt)[1] <- 'fustat'#给生存状态和生存时间改个代码
colnames(rt)[2] <- 'futime'

#rt=read.table("data.exp.txt",header=T,sep="\t",row.names=1)           
rt$futime=rt$futime/365   


set.seed(12)   #设置随机数种子，如果不设置，再来进行计算结果会不一样，同一个数字，可以保证计算结果相同
x=as.matrix(rt[,c(3:ncol(rt))])#基因的表达量
y=data.matrix(Surv(rt$futime,rt$fustat))#生存状况
fit=glmnet(x, y, family = "cox", maxit = 2000)
plot(fit, xvar = "lambda", label = TRUE)

cvfit = cv.glmnet(x, y, family="cox", maxit = 2000)
plot(cvfit)
#其中两条虚线分别指示了两个特殊的λ值
dev.off()
###4. 输出预测模型的相关系数与riskScore
###4.1 输出相关系数
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
geneCoef   #查看模型的相关系数#筛选出来的基因就是最佳的，计算风险评分：基因乘以对应系数再相加
write.table(geneCoef,file="lasso_gene.txt",sep="\t",quote=F,col.names=T)

###4.2 计算riskScore
FinalGeneExp = rt[,lassoGene]
myFun = function(x){crossprod(as.numeric(x),actCoef)}
riskScore = apply(FinalGeneExp,1,myFun)
outCol = c("futime", "fustat", lassoGene)
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))#进行判断，如果得出的分数大于中位数就是高组
dat = cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)#这个dat是独立的一个示例数据

###5. 绘制散点分布图
#install.packages("ggpubr")
library(ggpubr)  
p <- ggboxplot(dat, x = "fustat", y = "riskScore",#x代表生存状态
               color = "fustat", palette = "jco",
               add = "jitter")
p <- p + stat_compare_means()   #  Add p-value
p   #得出预测结果

###6. 判断预测结果的准确性
#有分组信息、风险评分了，之前讲的ROC是用基因表达来预测患者的生存，这个是用算出来的风险评分来预测
#install.packages("ROCR")
library(ROCR)   #使用ROCR包绘制预测模型的ROC曲线
library(glmnet)
library(caret)

pred <- prediction(dat$riskScore, dat$fustat)
perf <- performance(pred,"tpr","fpr")
AUC <- performance(pred,"auc")   #计算AUC
plot(perf,colorize=FALSE, col="red", print.auc =TRUE) #绘制ROC曲线
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
dev.off()

#画风险分布图
#y=生存时间
rt <- dat#这里进行赋予只是因为下面的代码都是rt
color=as.vector(rt$fustat)
color[color==1]="indianred1"
color[color==0]="lightseagreen"
plot(rt$futime, pch=19,
     xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("indianred1","lightseagreen"),cex=1.2)
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
abline(v=lowLength,lty=2)
dev.off()

#y=riskscore
rt <- rt[order(rt[,7]),]#通过rt数据框的第25列（即风险评分）进行排序（此处默认升序，如果想降序排列，在逗号后面加decreasing=t）
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
plot(line, type="p", pch=20,
     xlab="Patients (increasing risk socre)", ylab="Risk score",
     col=c(rep("lightseagreen",lowLength),rep("indianred1",highLength)) )
abline(h=median(rt$riskScore),v=lowLength,lty=2)#风险评分的中位数、低组长度
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("indianred1","lightseagreen"),cex=1.2)
dev.off()


#多因素coxgene可视化热图
#读取risk_train,在读进来之前用excel给高低风险排序
counts <- risk_train
counts <- counts[,-1]
counts <- counts[,-1]
counts <- counts[,-1]
counts1 <- counts[,-ncol(counts)] 
counts1 <- counts1[,-ncol(counts1)]
counts1 <- t(counts1)


library(pheatmap) # 加载pheatmap这个R包
#group_list=factor(ifelse(substr(counts$risk) == "high","high","low"),levels = c("low","high"))
group_list=counts$risk
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=rownames(counts)
#mydata1[,c(3:7)] <- as.numeric(unlist(mydata1[,c(3:7)]))
pheatmap(counts1,
         annotation_col=annotation_col,#需要的文件有表达谱exp-diff,分组信息group-list
         scale = "row",
         show_rownames = T,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         fontsize = 10,
         fontsize_row=10,
         fontsize_col=3)
dev.off()
