library(pheatmap) # 加载包
library(ggplot2) # 加载包

comgene <- intersect(rownames(risk_all),rownames(TCGA_KIRC_phenotype_sur))#取交集
table(substr(comgene,14,16))
risk_all <- risk_all[comgene,]
TCGA_KIRC_phenotype_sur <- TCGA_KIRC_phenotype_sur[comgene,]
#risk_all <- risk_all[,-5]#删掉不用的列
#TCGA_LUSC_phenotype_sur <- TCGA_LUSC_phenotype_sur[,-7]
data <- cbind(risk_all,TCGA_KIRC_phenotype_sur)
write.csv(data,file = "data.csv")
#把data的高低风险进行排序
data <- read.csv("C:/Users/caichenxi/Desktop/实践2肾癌 - 再次尝试/16.临床分期可视化热图/data.csv", row.names=1)


annotation_col = data.frame(risk = factor(data$risk),
                            age=factor(data$age),
                            gender=factor(data$gender),
                            stage=factor(data$stage))
data1 <- data[,(4:11)]
data1 <- t(data1)

rownames(annotation_col)
# [1] "1" "2" "3" "4" "5" "6" "7" "8" "9"
rownames(data)
# [1] "Q1" "Q2" "Q3" "F1" "F2" "F3" "T1" "T2" "T3"
rownames(annotation_col) <- rownames(data)
head(annotation_col)
#    Deal_with  Day
# Q1      Salt Day1
# Q2   Drought Day2
# Q3      Heat Day3
# F1      Salt Day1
# F2   Drought Day2
# F3      Heat Day3

p <- pheatmap(data1,scale="row",
              annotation_col = annotation_col,
              border="white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #去掉横、纵坐标id
              show_colnames = F,
              legend = T, # 添加图例
              legend_breaks=c(-2,0,2)) # 设置图例范围

pheatmap(data1,
         annotation_col=annotation_col,
         scale = "row",
         cluster_cols = F, # 去掉横向、纵向聚类#如果是T就不能高低风险排序
         cluster_rows = F,
         show_rownames = T,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize = 10,
         fontsize_row=10,
         fontsize_col=10)#画热图的代码
dev.off()
