
####肿瘤突变负荷TMB####
#TCGA突变数据下载
#网址：https://portal.gdc.cancer.gov/
library(TCGAbiolinks)#这个包要用BioManager装
#安装maftools
if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")
library(maftools)
library(tidyverse)

library(readxl)
library(readr)
mut2 <- load("kirc_maf.Rdata")
#突变数据下载后，整理成数据框
#mut2 <- read.maf("TCGA.LIHC.varscan.40fe9c1b-19d0-45cf-898a-f7b0cbad783e.DR-10.0.somatic.maf")
mut2 <- KIRC.maf
a <- mut2@data %>% 
  .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>% 
  as.data.frame() %>% 
  mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,12))

gene <- as.character(unique(a$Hugo_Symbol))#把a中第一列的基因提取出来，进行去重复，并变为字符串的形式
sample <- as.character(unique(a$Tumor_Sample_Barcode))

mat <- as.data.frame(matrix("",length(gene),length(sample),#跑出来的mat表里有的突变，在下面mat_0_1里就表示为1
                            dimnames = list(gene,sample)))
mat_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                dimnames = list(gene,sample)))

for (i in 1:nrow(a)){
  mat[as.character(a[i,1]),as.character(a[i,3])] <- as.character(a[i,2])
}
for (i in 1:nrow(a)){
  mat_0_1[as.character(a[i,1]),as.character(a[i,3])] <- 1
}

gene_count <- data.frame(gene=rownames(mat_0_1),#列表，看每种基因有多少次突变
                         count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))
gene_top <- gene_count$gene[1:20] # 修改数字，代表TOP多少
save(mat,mat_0_1,file = "TMB-KIRC.rda")

#以下绘图代码来自解螺旋阿琛老师
oncoplot(maf = mut2,#这里用的数据是mut2，是刚才在上面读取的数据
         top = 30,   #显示前30个的突变基因信息
         fontSize = 0.6,   #设置字体大小
         showTumorSampleBarcodes = F)   #不显示病人信息
dev.off()

####计算TMB####
maf = tmb(maf = mut2,
          captureSize = 50,
          logScale = TRUE)   
maf$sample <- substr(maf$Tumor_Sample_Barcode,1,16)
maf$sample <- gsub("-",".",maf$sample)
rownames(maf) <- maf$sample
#手动更改后导入
write.csv(maf, file = "maf.csv")
group <- risk_list
rownames(group) <- substring(rownames(group),1,12)
write.csv(group,file = "risk_list1.csv")
com <- intersect(rownames(maf),rownames(group))
maf <- maf[com,]
group$x <- as.character(group$V2)
group <- group %>% t() %>% as.data.frame()
group <- group[-1,]
group <- group[,com]
group <- group %>% t() %>% as.data.frame()
identical(rownames(group),rownames(maf))
tmb <- cbind(maf,group)
write.csv(tmb, file = "tmb.csv")
