## 组间相关性分析
## 皮尔森相关性（pearson)衡量两个连续变量的线性相关关系；
## 斯皮尔曼相关性(spearman)衡量两个变量之间的单调关系，两个变量同时变化，但不一定是线性关系；


setwd("C:/Users/24925/Desktop/Rdata")
library(psych)
library(ggplot2)
library(pheatmap)

comgene <- intersect(rownames(KIRC_fpkm_all),rownames(ICDgene))#intersect取交集
relate_gene_exp <- KIRC_fpkm_all[comgene,]
comgene <- intersect(colnames(TCGA_KIRC_lncRNA),colnames(relate_gene_exp))
TCGA_KIRC_lncRNA <- TCGA_KIRC_lncRNA[,comgene]
relate_gene_exp <- relate_gene_exp[,comgene]
write.csv(relate_gene_exp, "ICDgene_exp.csv",row.names = T)
write.csv(TCGA_KIRC_lncRNA, "TCGA_KIRC_lncRNA_all.csv",row.names = T)

#TCGA_COADREAD_lncRNA_a <- TCGA_COADREAD_lncRNA[(1:7415),]
#TCGA_COADREAD_lncRNA_b <- TCGA_COADREAD_lncRNA[(7416:14831),]
#write.csv(TCGA_COADREAD_lncRNA_a, "TCGA_COADREAD_lncRNA_a.csv",row.names = T)
#write.csv(TCGA_COADREAD_lncRNA_b, "TCGA_COADREAD_lncRNA_b.csv",row.names = T)

#diff_exp_cor <- rbind(p_list1,p_list2)



Diffsing_corr <- diff_exp_cor[(diff_exp_cor$p_value< 0.05 & abs(diff_exp_cor$cor) > 0.3),]
library(limma)
Diffsing_corr = as.data.frame(avereps(Diffsing_corr[,-2],ID = Diffsing_corr$Var1) )
write.csv(Diffsing_corr,"Diffsing.Corr.output.csv", row.names = F)







# 导入数据格式01
# 直接可以只使用，无需进行转换
df01 <- read.table('m7g.txt', header = T, sep = '\t')
df02 <- read.table('lncrna.txt',header = T,sep = '\t')

## 导入数据格式02，需要进行装换
#df03 <- read.table("lncRNA_inputdata_格式2.txt", row.names = 1, header = T,sep = '\t')
#df04 <- read.table("mRNA_inputdata_格式2.txt", row.names = 1, header = T,sep = '\t')
# 转换数据格式
## 样本名称保持一致
df03 <- TCGA_LUSC_lncRNA_all
df04 <- relate_gene_exp
#df03 <- df03[,colnames(df03)]

## 如果在准备数据时，样本名称不一致可以按以下操作
 #df03 <- df03[,colnames(df04)]

df03 <- data.frame(t(df03))
df04 <- df04[,colnames(df04)]
df04 <- data.frame(t(df04))
df03[1:4,1:4]
df04[1:4,1:4]

## --------
## 过滤 表达量都为0的
df03 <- df03[,apply(df03,2,sum) > 0]
df04 <- df04[,apply(df04,2,sum) > 0]
## ------------------
## 使用 corr.test()进行计算
corr <- corr.test(df03, df04, use = 'pairwise', 
                 method = 'pearson', adjust = 'BH', ## fdr
                 alpha = 0.05)
save(corr,file = "corr.rda")
corr$p.adj
corr$r
##-----------------
## 导出相关的Corr值和P值
###-----------方法一 
pdat <- corr$r  ## 相关性
corr.data <- data.frame(x = rep(colnames(pdat), each = 49),
                        y = rep(rownames(pdat),each=14357),
                        Cor = as.vector(corr$r),
                        P = as.vector(corr$p.adj))
head(corr.data)
write.csv(corr.data, "Corr.output.csv",row.names = F)
## 筛选出符合自己要求的基因对
Diffsing_corr <- diff_exp_cor[(diff_exp_cor$p_value< 0.001 & abs(diff_exp_cor$cor) > 0.5),]
library(limma)
Diffsing_corr = as.data.frame(avereps(Diffsing_corr[,-2],ID = Diffsing_corr$Var1) )
write.csv(Diffsing_corr,"Diffsing.Corr.output.csv", row.names = F)


## ----     绘制相关性图
## 热图 一
library(pheatmap)
pdat <- corr$r  ## 相关性
pdat2 <- matrix(ifelse(abs(corr$r) <= 0.3, 0, pdat), nrow(pdat)) ## corr小于0.3的变为0,nrow(pdat):规定向量的行
pdat2 <- matrix(ifelse(corr$p.adj >= 0.05, 0, pdat2), nrow(pdat))  ## P值小于大于0.05，变为0

colnames(pdat2) <- colnames(pdat)
rownames(pdat2) <- rownames(pdat)

bk <- seq(-1, 1, by = 0.1) ## 刻度的分割
bk
length(bk)

# 颜色的设置
mycol <- c(colorRampPalette(c('blue', 'white'))(9), 
           colorRampPalette(c('white', 'white'))(3),
           colorRampPalette(c('white', 'red'))(9))
pdf("one_相关热图.pdf", width = 6, height = 6)
pheatmap(pdat2, cellwidth = 1, cellheight = 1,
         cluster_cols = F,
         color = mycol,
         legend_breaks = seq(-1,1,by = 0.2),
         breaks = bk,
         angle_col = 45)
dev.off()

### 气泡图
## 导出数据
pdat3 <- data.frame(x = rep(colnames(pdat2), each = 49), ## y纵坐标复制 23遍
                    y = rep(rownames(pdat2), 14357),        ## x横坐标复制13遍
                    Cor = as.vector(corr$r),               ## 赋值corr值
                    P = as.vector(
                      ifelse(corr$p.adj >= 0.05, 1, corr$p.adj)
                    ))
head(pdat3)
library(ggplot2)
pdf("one_气泡图.pdf",width = 6, height = 6)
ggplot(pdat3, aes(x,y, colour = Cor, size = -log10(P))) +
  geom_point(aes(alpha = -log10(P))) + 
  scale_color_gradient2(low = 'blue',
                        mid = 'white',
                        high = 'red',
                        midpoint = 0) +
  scale_size_area() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  ## 更改横纵坐标轴中的字体颜色、大小
  theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(colour = "black", size = 10))
dev.off()
