#设置工作位置
setwd("DataFiles")

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(limma)
library(ggplot2)
#设置训练集
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='Training Data'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

#读入表达量数据，行名为基因，列名为样本不能有重复样本名
lncrnaexp<- read.table("hualiao.txt",header=T,sep="\t")
lncrnaexp<-`TCGA_KIRC_fpkm_log2+1`
#rownames(lncrnaexp)=lncrnaexp[,1]
#lncrnaexp <- lncrnaexp[,-1]
#lncrnaexp = t(lncrnaexp)
tumor<- colnames(lncrnaexp)[substr(colnames(lncrnaexp),14,15) == "01"]#提取肿瘤样本
lncrnaexp <- lncrnaexp[,tumor]
lncrnaexp = as.matrix(lncrnaexp)

#药物预测
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = lncrnaexp,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 1, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )


#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#风险文件和药物敏感性合并
sameSample=intersect(row.names(risk_all), rownames(DrugPredictions))
risk=risk_all[sameSample, "risk",drop=F]
DrugPredictions=DrugPredictions[sameSample,]
rt=cbind(risk, DrugPredictions)

#设置比较组
rt$risk=factor(rt$risk, levels=c("low", "high"))
type=levels(factor(rt[,"risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


##test
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(rt)]
for(i in colnames(rt)[1:(ncol(rt)-1)]){boxplot=ggboxplot(rt, x="risk", y=i, fill="risk",
                                                         xlab="Risk",
                                                         ylab=i,
                                                         legend.title="Risk",
                                                         palette=c("green", "red"))+ 
  stat_compare_means(comparisons=my_comparisons)
pdf(file=paste0("boxplot.", i, ".pdf"), width=6, height=5.5)
print(boxplot)
dev.off()}





#绘制箱线图
boxplot=ggboxplot(rt, x="risk", y="Dactolisib_1057", fill="risk",#此处更改药物名称
                  xlab="Risk",
                  ylab= "Dactolisib_1057 senstivity (IC50)",
                  legend.title="Risk",
                  palette=c("green", "red")
)+ 
  stat_compare_means(comparisons=my_comparisons)
pdf(file="Dactolisib_1057 sentivity.pdf", width=5, height=4.5)
print(boxplot)
dev.off()
