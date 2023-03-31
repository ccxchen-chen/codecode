

#数据整理，差异表达的基因（pearson结果）#
KIRC_fpkm_all <- TCGA_KIRC_fpkmdata_log2+1
comgene <- intersect(rownames(Diffsing.Corr.output),rownames(KIRC_fpkm_all))#取交集
table(substr(comgene,14,16))
diff_gene_exp <- KIRC_fpkm_all[comgene,]
write.csv(diff_gene_exp,file = "diff_gene_exp_fpkm.csv")

#整理全部基因，注意要用risk_all取交集，去掉正常人
KIRC_fpkm_all <- as.data.frame(t(KIRC_fpkm_all))
comgene <- intersect(rownames(risk_all),rownames(KIRC_fpkm_all))#取交集
table(substr(comgene,14,15))
KIRC_fpkm_all <- KIRC_fpkm_all[comgene,]
KIRC_fpkm_all_risk=cbind(risk_all$risk,KIRC_fpkm_all)
colnames(KIRC_fpkm_all_risk)[1] <- 'risk'
write.csv(KIRC_fpkm_all_risk,file = "KIRC_fpkm_all_risk.csv")

#整理lncRNA
KIRC_fpkm_all <- KIRC_lncRNA_all
KIRC_fpkm_all <- as.data.frame(t(KIRC_fpkm_all))
comgene <- intersect(rownames(risk_all),rownames(KIRC_fpkm_all))#取交集
table(substr(comgene,14,15))
KIRC_fpkm_all <- KIRC_fpkm_all[comgene,]
KIRC_lncRNA_all_risk=cbind(risk_all$risk,KIRC_fpkm_all)#在跑之前注意改变前面的文件名
colnames(KIRC_lncRNA_all_risk)[1] <- 'risk'
write.csv(KIRC_lncRNA_all_risk,file = "KIRC_lncRNA_all_risk.csv")

#整理差异表达基因
diff_gene_exp<- as.data.frame(t(diff_gene_exp))
KIRC_fpkm_all <- diff_gene_exp
comgene <- intersect(rownames(risk_all),rownames(KIRC_fpkm_all))#取交集
table(substr(comgene,14,15))
KIRC_fpkm_all <- KIRC_fpkm_all[comgene,]
diff_gene_exp_risk=cbind(risk_all$risk,KIRC_fpkm_all)#在跑之前注意改变前面的文件名
colnames(diff_gene_exp_risk)[1] <- 'risk'
write.csv(diff_gene_exp_risk,file = "diff_gene_exp_risk.csv")

#mcoxgene数据整理
risk_all <- risk_all[,-1]
risk_all <- risk_all[,-1]
risk_all <- risk_all[,-1]
risk_all <- risk_all[,-9]
# 得到 dataframe 的列名数组
cols <- colnames(risk_all)
# 根据需要，生成新的列名顺序，例如，把倒数第一列插入到正数第二列之前，假设目前的列名顺序是
# A B C D E F G H
# 操作以后会变成
# A H B C D E F G
new_cols <- c(cols[length(cols)], cols[1:(length(cols) - 1)])
# 然后将 dataframe 按照新的列名顺序排列
risk_all <- risk_all[, new_cols]
write.csv(risk_all,file = "risk_all_PCA.csv")
#install.packages("scatterplot3d")


library(scatterplot3d)      
   
#outFile="PCA.pdf"          
#setwd("D:\\biowolf\\bioR\\46.3dPCA")    

rt=KIRC_fpkm_all_risk
data=rt[,c(2:ncol(rt))]
Type=rt[,1]  #提取分组信息
var=colnames(rt)[1]

#颜色
group=levels(factor(Type))
bioCol=c("red","blue","green","yellow")
col= bioCol[match(Type,group)]

#去掉表达全部为0的列#从网上找的
del <- c()#定义向量来存储值
for (i in seq(1, ncol(data))) {
  if(sum(data[,i])==0){
    #print(i)
    del <- append(del, -i)
  }
}
data <- data[,del]
#PCA分析
data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)

#绘制
#pdf(file=outFile, height=5, width=6)
par(oma=c(0.5,0.5,0.5,0.5))
s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=col)
legend("top", legend =group,pch = 16, inset = -0.2, box.col="white",xpd = TRUE, horiz = TRUE,col=bioCol[1:length(group)])
dev.off()



