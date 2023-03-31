#免疫浸润分析
#数据要求：使用FPKM和TPM或DESEq2标准化后的矩阵，行名为基因名，列为样本名
#需要的数据为肿瘤mRNA数据，由mRNA算出免疫细胞浸润
#setwd("C:/Users/24925/Desktop/Rdata")

remove(list = ls()) #一键清空
#加载包
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)

source('Cibersort.R')

# 设置分析依赖的基础表达文件
# 每类免疫细胞的标志性基因及其表达
# 基因名字为Gene symbol
LM22.file <- "LM22.txt"
write.table(TCGA_LUSC_mRNA, file = "TCGA_LUSC_mRNA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
TCGA_LUSC_mRNA <- read.delim("C:/Users/caichenxi/Desktop/test/15.cibersort/TCGA_LUSC_mRNA.txt", row.names=1)
write.table(TCGA_LUSC_mRNA_all, file = "TCGA_LUSC_mRNA_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
# 1. Cibersort

TCGA_exp.file <- "TCGA_KIRC_mRNA_all.txt"

TCGA_TME.results <- CIBERSORT(LM22.file ,TCGA_exp.file, perm = 50, QN = F)  #不可以用csv文件打开
# perm置换次数=1000
# QN如果是芯片设置为T，如果是测序就设置为F

write.csv(TCGA_TME.results, "TCGA_CIBERSORT_Results_fromRcode.csv")


# 2. 分组信息
# TCGA的数据可以从名字获取
TCGA_TME.results <- as.data.frame(TCGA_TME.results)#重点！！！
comgene <- intersect(rownames(TCGA_TME.results),rownames(risk_all))#intersect取交集

TCGA_TME.results1 <- TCGA_TME.results[comgene,]
risk_all <- risk_all[comgene,]#记得要给两个都改成取交集之后的
identical(rownames(TCGA_TME.results1),rownames(risk_all))
TCGA_TME.results1 <- cbind(TCGA_TME.results1,risk_all$risk)
colnames(TCGA_TME.results1)[26] <- 'risk'


group_list <- factor(TCGA_TME.results1$risk, levels = c("low","high"))



table(group_list) # Normal    Tumor 

## group_list
## Nontumor    Tumor 
             

## 3. 绘图（箱线图）
# 3.1 数据粗处理
TME_data <- as.data.frame(TCGA_TME.results1[,1:22])

TME_data$group <- group_list
TME_data$sample <- row.names(TME_data)

# 3.2 融合数据
TME_New = melt(TME_data)

## Using group, sample as id variables

colnames(TME_New)=c("Group","Sample","Celltype","Composition")  #设置行名
head(TME_New)

##      Group          Sample      Celltype Composition
## 1    Tumor TCGA.CV.6943.01 B cells naive 0.007651678
## 2    Tumor TCGA.CV.6959.01 B cells naive 0.019549031
## 3 Nontumor TCGA.CV.7438.11 B cells naive 0.025349204
## 4 Nontumor TCGA.CV.7242.11 B cells naive 0.032583659
## 5    Tumor TCGA.CV.7432.01 B cells naive 0.000000000
## 6 Nontumor TCGA.CV.6939.11 B cells naive 0.074282293

# 3.3 按免疫细胞占比中位数排序绘图（可选）
#plot_order = TME_New[TME_New$Group=="Tumor",] %>% 
#  group_by(Celltype) %>% 
#  summarise(m = median(Composition)) %>% 
#  arrange(desc(m)) %>% 
#  pull(Celltype)

## `summarise()` ungrouping output (override with `.groups` argument)

#TME_New$Celltype = factor(TME_New$Celltype,levels = plot_order)


# 3.3 出图
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }

box_TME <- ggplot(TME_New, aes(x = Celltype, y = Composition))+ 
  labs(y="Cell composition",x= NULL,title = "TME Cell composition")+  
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c( "#1CB4B8","#EB7369"))+
  
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)

box_TME;ggsave("TCGA_KIRC_TME.pdf",box_TME,height=15,width=20,unit="cm")
