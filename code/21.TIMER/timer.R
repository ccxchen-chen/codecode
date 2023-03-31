library(tidyverse)
library(pheatmap) # 加载包
library(ggplot2) # 加载包
comname <- intersect(colnames(TCGA_KIRC_mRNA_all),rownames(risk_all)) 
TCGA_KIRC_mRNA_01 <- TCGA_KIRC_mRNA_all[,comname]
write.csv(TCGA_KIRC_mRNA_01,file = "TCGA_KIRC_mRNA_01.csv")

TCGA_KIRC_mRNA_01a <- TCGA_KIRC_mRNA_01[,(1:258)]
write.csv(TCGA_KIRC_mRNA_01a,file = "TCGA_KIRC_mRNA_01a.csv")

TCGA_KIRC_mRNA_01b <- TCGA_KIRC_mRNA_01[,(259:526)]
write.csv(TCGA_KIRC_mRNA_01b,file = "TCGA_KIRC_mRNA_01b.csv")

#TCGA_KIRC_mRNA_01c <- TCGA_KIRC_mRNA_01[,(321:493)]
#write.csv(TCGA_LUSC_mRNA_01c,file = "TCGA_LUSC_mRNA_01c.csv")


risk_all <- risk_all[comname,]
TCGA_KIRC_mRNA_01 <- TCGA_KIRC_mRNA_01 %>% t() %>% as.data.frame()
TCGA_KIRC_mRNA_01$risk <- as.factor(risk_all$risk)
TCGA_KIRC_mRNA_01$riskScore <- as.factor(risk_all$riskScore)
write.csv(TCGA_KIRC_mRNA_01,file = "TCGA_KIRC_mRNA_01_risk.csv")

#读取timer网站下载文件
timer <- cbind(estimation_01matrix,estimation_02matrix)
timer <- timer %>% t() %>% as.data.frame()

comname <- intersect(rownames(TCGA_KIRC_mRNA_01_risk),rownames(timer)) 
TCGA_KIRC_mRNA_01_risk <- TCGA_KIRC_mRNA_01_risk[comname,]
timer <- timer[comname,]
timer$risk <- as.factor(TCGA_KIRC_mRNA_01_risk$risk)
timer$riskScore <- as.factor(TCGA_KIRC_mRNA_01_risk$riskScore)
write.csv(timer,file = "timer_risk.csv")

data <- timer_risk_less
#annotation_col = data.frame(risk = factor(data$risk),
                            #riskScore=factor(data$riskScore))
annotation_col = data.frame(risk = factor(data$risk))
data1 <- data[,(1:98)]
data1 <- t(data1)

rownames(annotation_col)
# [1] "1" "2" "3" "4" "5" "6" "7" "8" "9"
rownames(data)
# [1] "Q1" "Q2" "Q3" "F1" "F2" "F3" "T1" "T2" "T3"
rownames(annotation_col) <- rownames(data)
head(annotation_col)
ann_colors = list(risk = c("low"=  "blue", "high"=  "red")) #定义分组颜色
#    Deal_with  Day
# Q1      Salt Day1
# Q2   Drought Day2
# Q3      Heat Day3
# F1      Salt Day1
# F2   Drought Day2
# F3      Heat Day3

annotation_row = data.frame(methods = factor(grouplist_less$methods))
rownames(annotation_row) <- rownames(data1)
#head(annotation_row)
#       GeneClass
# Gene1      WRKY
# Gene2       AP2
# Gene3     YABBY
# Gene4      WRKY
# Gene5       AP2
# Gene6     YABBY

pheatmap(data, annotation_row =annotation_row)



#ann_colors = list(risk = c("low"=  "#f3bf88", "high"=  "#f7b977")) #定义分组颜色
pheatmap(data1,
         annotation_col=annotation_col,
         annotation_row =annotation_row,
         scale = "row",
         cluster_cols = F, # 去掉横向、纵向聚类#如果是T就不能高低风险排序
         cluster_rows = F,
         show_rownames = T,
         show_colnames =F,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         fontsize = 5,
         fontsize_row=4,
         fontsize_col=5)#画热图的代码
dev.off()

