####cox回归分析####
#设置工作目录
setwd("cox")
#安装加载R包
install.packages("survival")
install.packages("forestplot")
library(survival)
library(forestplot)
library(tidyverse)
#下载生存信息
#xena官网：https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Liver%20Cancer%20(LIHC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#读取生存信息tsv文件
surv = read.table(file = 'TCGA-KIRC.survival.tsv', sep = '\t', header = TRUE) 

#整理生存信息数据
surv$sample <- substring(surv$sample,1,15) %>% gsub("-",".",.)
surv <- surv[!duplicated(surv$sample),]   #去重复

rownames(surv) <- surv$sample

surv <- surv[,-1]
surv <- surv[,-2]
#读取表达数据

expr <- read.table("TCGA_KIRC_fpkm.csv",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
expr <- TCGA_KIRC_fpkm
rownames(expr) <- expr[,1]
expr[1:4,1:4]
###
tumor <- colnames(expr)[substr(colnames(expr),14,15) == "01"]
counts_01A <- expr[,tumor]
normal <- colnames(expr)[substr(colnames(expr),14,15) == "11"]
counts_11A <- expr[,normal]
TCGA_KIRC_fpkm <- cbind(counts_11A,counts_01A)
write.csv(TCGA_KIRC_fpkm,file = "TCGA_KIRC_fpkm_TN.csv")
expr <- TCGA_KIRC_fpkm

#expr <- ceiling(2^(expr)-1)
#write.csv(expr,file = "TCGA_COADREAD_fpkm_nlog.csv")
comgene <- intersect(colnames(expr),rownames(surv))#取交集

table(substr(comgene,14,15))#分组计数
expr <- expr[,comgene]#改列名
surv <- surv[comgene,]

#表达数据整理完毕
#读取tcga差异分析结果（DEG）
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 3, padj < 0.05)#此处认为超过两倍为有差异性
#整合
deg_expr <- expr[rownames(res_deseq2),] %>% t() %>% as.data.frame()
surv.expr <- cbind(surv,deg_expr)#cbind表示按照列进行拼接；rbind表示按照行进行拼接
write.csv(surv.expr,file = "ucox_surv.expr_all.csv")
write.table(surv.expr, file = "ucox_surv.expr_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#去除正常患者
tumor <- rownames(ucox_surv.expr_all)[substr(rownames(ucox_surv.expr_all),14,15) == "01"]
ucox_surv.expr_01 <- ucox_surv.expr_all[tumor,]
write.csv(ucox_surv.expr_01,file = "ucox_surv.expr_01.csv")
write.table(ucox_surv.expr_01, file = "ucox_surv.expr_01.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



####分组####
#(脚本)
install.packages("caret")
setwd("C:/Users/caichenxi/Desktop/test/4.单因素cox")
rt=read.table("ucox_surv.expr_01.txt",sep = "\t",header=T,check.names = F)

rownames(rt) <- rt[,1]  
rt = rt[,-1]

library(caret)
inTrain <-createDataPartition(y=rt[,3],p=0.5,list=F)#p=0.5表示训练组有50%
train <- rt[inTrain,]
test <- rt[-inTrain,]
write.table(train,file="ucox_surv.expr_train.txt",sep="\t",quote=F,row.names=T)
write.table(test,file="ucox_surv.expr_test.txt",sep="\t",quote=F,row.names=T)

write.csv(train, file = "ucox_surv.expr_train.csv")
write.csv(test, file = "ucox_surv.expr_test.csv")

#Cox分析
surv.expr1=read.table("ucox_surv.expr_train.txt",sep = "\t",header=T,check.names = F)
#surv.expr1 <- ucox_surv.expr_train
rownames(surv.expr1) <- surv.expr1[,1]  
#surv.expr1 = surv.expr1[,-1]

Coxoutput <- NULL 
for(i in 3:ncol(surv.expr1)){
  g <- colnames(surv.expr1)[i]
  cox <- coxph(Surv(OS.time,OS) ~ surv.expr1[,i], data = surv.expr1)# 单变量cox模型
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}


write.table(Coxoutput, file = "cox results_train.txt",sep = "\t",row.names = F,col.names = T,quote = F)
###筛选top基因
pcutoff <- 0.05#p值的截断值
topgene <- Coxoutput[which(Coxoutput$pvalue < pcutoff),]# 取出p值小于阈值的基因
write.csv(topgene, file = "topgene.csv")
topgene <- topgene[order(topgene[,4]),]#通过rt数据框的第25列（即风险评分）进行排序（此处默认升序，如果想降序排列，在逗号后面加decreasing=t）
write.csv(topgene,file = "cox_topgene_01.csv")
topgene <- topgene[1:10,]#提取前十个


#3. 绘制森林图
##3.1 输入表格的制作
tabletext <- cbind(c("Gene",topgene$gene),
                   c("HR",format(round(as.numeric(topgene$HR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(topgene$lower),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(topgene$upper),3),nsmall = 3)),
                   c("pvalue",format(round(as.numeric(topgene$p),3),nsmall = 3)))
##3.2 绘制森林图
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(topgene$HR)),
           lower=c(NA,as.numeric(topgene$lower)), 
           upper=c(NA,as.numeric(topgene$upper)),
           graph.pos=5,# 图在表中的列位置
           graphwidth = unit(.25,"npc"),# 图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",# box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),# box颜色
           
           boxsize=0.4,# box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=T,# 显示区间
           zero=1,# zero线横坐标
           lwd.zero=1.5,# zero线宽
           xticks = c(0.5,1,1.5),# 横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2),# 各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"), # 在第一行上面画黑色实线
                           "2" = gpar(lwd=1.5, col="black"), # 在第一行标题行下画黑色实线
                           "12" = gpar(lwd=2, col="black")), # 在最后一行上画黑色实线
           lineheight = unit(.9,"cm"),# 固定行高
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
dev.off()


#cox分析基因可视化热图
topgene <- cox_topgene_train
KIRC_fpkm_all <- TCGA_KIRC_fpkm
rownames(topgene) <- topgene[,1]  
topgene = topgene[,-1]
comgene <- intersect(rownames(topgene),rownames(KIRC_fpkm_all))#取交集
cox_topgene_exp <- KIRC_fpkm_all[comgene,]
counts <- cox_topgene_exp
tumor <- colnames(counts)[substr(colnames(counts),14,15) == "01"]
counts_01A <- counts[,tumor]
normal <- colnames(counts)[substr(colnames(counts),14,15) == "11"]
counts_11A <- counts[,normal]
counts <- cbind(counts_11A,counts_01A)


library(pheatmap) # 加载pheatmap这个R包
group_list=factor(ifelse(substr(colnames(counts),14,15) == "01","T","N"),levels = c("N","T"))

annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(counts)
pheatmap(counts,
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

