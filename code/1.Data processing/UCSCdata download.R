setwd("C:/Users/24925/Desktop/Rdata")
#### 载入R包
library(openxlsx)
library(tidyverse)

library(limma)
library(readr)

#### 读取表达数据，并做粗处理。可以输入Count数据，也可以用FPKM数据

TCGA_rawdata <- read_tsv("TCGA-KIRC.htseq_counts.tsv")
TCGA_rawdata <- TCGA.KIRC.htseq_fpkm
dim(TCGA_rawdata) 

###读取基因标注文件
probeMap <- read.table("gencode.v22.annotation.gene.probeMap",sep = "\t" , header = T)
probeMap[1:4,1:4]

###转换ID
TCGA_gset <- TCGA_rawdata %>%
  inner_join(probeMap, by = c("Ensembl_ID" = "id")) %>%
  select(gene, starts_with("TCGA") )
TCGA_gset[1:4,1:4]





###利用limma包对重复的基因取平均值合并
TCGA_gset = as.data.frame(avereps(TCGA_gset[,-1],ID = TCGA_gset$gene) )

###列名重命名：取列名的1-15位
colnames(TCGA_gset) <- substring(colnames(TCGA_gset),1,15) %>% gsub("-",".",.)
write.csv(TCGA_gset,"TCGA_KIRC_counts_log2+1.csv")
TCGA_gset[1:4,1:4]

###根据表达矩阵病人的ID进行分组
TCGA_group_list <- ifelse(as.numeric(substring(colnames(TCGA_gset),14,15)) < 10,
                          "Tumor","Normal") %>% 
  factor(.,levels = c("Normal","Tumor"))
table(TCGA_group_list)    

### 注释mRNA，lncRNA和miRNA
mRNA_info <- read.xlsx("Gene_info.xlsx",sheet = "mRNA_info")
lncRNA_info <- read.xlsx("Gene_info.xlsx",sheet = "lncRNA_info")
#miRNA_info <- read.xlsx("Gene_info.xlsx",sheet = "miRNA_info")

###根据基因的注释信息，提取对应的表达矩阵
mRNA_gset <- TCGA_gset[rownames(TCGA_gset) %in% mRNA_info$gene_name,]
dim(mRNA_gset) 

write.csv(mRNA_gset,"TCGA_KIRC_mRNA_all.csv",quote = F,row.names = T)

lncRNA_gset <- TCGA_gset[rownames(TCGA_gset) %in% lncRNA_info$gene_name,]
dim(lncRNA_gset) 

write.csv(lncRNA_gset,"TCGA_KIRC_lncRNA_all.csv",quote = F,row.names = T)

miRNA_gset <- TCGA_gset[rownames(TCGA_gset) %in% miRNA_info$gene_name,]
dim(miRNA_gset) 

write.csv(miRNA_gset,"TCGA_LUAD_miRNA.csv",quote = F,row.names = T)




#### 3. 读取表型数据，并做粗处理，用来做后续分析，比如生存分析
Phenodata <- read_tsv("TCGA-KIRC.GDC_phenotype.tsv")

Phenodata[1:4,1:4]

Phenodata$submitter_id.samples <- substring(Phenodata$submitter_id.samples,1,15) %>% 
  gsub("-",".",.)
Phenodata[1:4,1:4]

###survival数据
Sur_data <- read_tsv("TCGA-KIRC.survival.tsv")

Sur_data$sample <- substring(Sur_data$sample,1,15) %>% gsub("-",".",.)
Sur_data[1:4,1:4]

###结合，提取部分列
Phen_surv <- Phenodata %>%
  inner_join(Sur_data,by = c("submitter_id.samples" = "sample")) %>%
  select(submitter_id.samples,age_at_index.demographic,gender.demographic,race.demographic,
         pathologic_T,pathologic_N,pathologic_M,tumor_stage.diagnoses,
         vital_status.demographic,OS,OS.time)
head(Phen_surv)

###表型数据和表达数据匹配并排序,并加上分组
TCGA_gset <- t(TCGA_gset)
TCGA_gset <- as.data.frame(TCGA_gset)
write.csv(TCGA_gset,"TCGA_KIRC_gset_all.csv")

Phen_surv = Phen_surv[match(rownames(TCGA_gset),Phen_surv$submitter_id.samples),]
#write.csv(Phen_surv,"Phen_surv.csv")
identical(Phen_surv$submitter_id.samples,rownames(TCGA_gset))

#identical是false，就进行下列操作让他变true#取行名交集
comgene <- intersect(rownames(TCGA_gset),Phen_surv$submitter_id.samples)#intersect取交集
TCGA_gset <- TCGA_gset[comgene,]
#Phen_surv <- Phen_surv[comgene,]
class(TCGA_gset)#判断数据类型
class(comgene)
Phen_surv$submitter_id.samples <- comgene
a <- rownames(TCGA_gset)
b <- Phen_surv$submitter_id.samples
identical(a,b)#判断两个或多个元素是否完全相同，包括顺序

## [1] TRUE

Phen_surv$group <- TCGA_group_list
Phen_surv = dplyr::select(Phen_surv,submitter_id.samples,group,everything())
write.csv(Phen_surv,"TCGA_KIRC_phenotype_sur.csv")

head(Phen_surv)



###数据合并
comgene <- intersect(rownames(`TCGA_COAD_fpkm_log2+1`),row.names(`TCGA_READ_fpkm_log2+1`))#intersect取交集
`TCGA_COAD_fpkm_log2+1` <- `TCGA_COAD_fpkm_log2+1`[comgene,]
`TCGA_READ_fpkm_log2+1` <- `TCGA_READ_fpkm_log2+1`[comgene,]
TCGA_COADREAD_fpkm <- cbind(`TCGA_COAD_fpkm_log2+1`,`TCGA_READ_fpkm_log2+1`)
write.csv(TCGA_COADREAD_fpkm,"TCGA_COADREAD_fpkm.csv")

comgene <- intersect(rownames(`TCGA_COAD_counts_log2+1`),row.names(`TCGA_READ_counts_log2+1`))#intersect取交集
`TCGA_COAD_counts_log2+1` <- `TCGA_COAD_counts_log2+1`[comgene,]
`TCGA_READ_counts_log2+1` <- `TCGA_READ_counts_log2+1`[comgene,]
TCGA_COADREAD_counts <- cbind(`TCGA_COAD_counts_log2+1`,`TCGA_READ_counts_log2+1`)
write.csv(TCGA_COADREAD_counts,"TCGA_COADREAD_counts.csv")

comgene <- intersect(rownames(TCGA_COAD_lncRNA_all),row.names(TCGA_READ_lncRNA_all))#intersect取交集
TCGA_COAD_lncRNA_all <- TCGA_COAD_lncRNA_all[comgene,]
TCGA_READ_lncRNA_all <- TCGA_READ_lncRNA_all[comgene,]
TCGA_COADREAD_lncRNA <- cbind(TCGA_COAD_lncRNA_all,TCGA_READ_lncRNA_all)
write.csv(TCGA_COADREAD_lncRNA,"TCGA_COADREAD_lncRNA.csv")

comgene <- intersect(rownames(TCGA_COAD_mRNA_all),row.names(TCGA_READ_mRNA_all))#intersect取交集
TCGA_COAD_mRNA_all <- TCGA_COAD_mRNA_all[comgene,]
TCGA_READ_mRNA_all <- TCGA_READ_mRNA_all[comgene,]
TCGA_COADREAD_mRNA <- cbind(TCGA_COAD_mRNA_all,TCGA_READ_mRNA_all)
write.csv(TCGA_COADREAD_mRNA,"TCGA_COADREAD_mRNA.csv")














