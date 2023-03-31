

library(tidyverse)
surv = read.table(file = 'TCGA-KIRC.GDC_phenotype.tsv', sep = '\t', header = TRUE) 
surv <- TCGA.KIRC.GDC_phenotype
#整理生存信息数据
surv$submitter_id.samples <- substring(surv$submitter_id.samples,1,15) %>% gsub("-",".",.)
surv <- surv[!duplicated(surv$submitter_id.samples),]   #去重复

rownames(surv) <- surv$submitter_id.samples

comgene <- intersect(rownames(surv),rownames(ucox_surv.expr_01))#取交集

table(substr(comgene,14,15))#分组计数

surv <- surv[comgene,]

write.csv(surv,file = "TCGA-KIRC.GDC_phenotype526.csv")

comgene <- intersect(rownames(surv),rownames(ucox_surv.expr_test))#取交集

table(substr(comgene,14,15))#分组计数

surv <- surv[comgene,]
write.csv(surv,file = "TCGA-KIRC.GDC_phenotype_test.csv")
comgene <- intersect(rownames(surv),rownames(ucox_surv.expr_train))#取交集

table(substr(comgene,14,15))#分组计数

surv <- surv[comgene,]
write.csv(surv,file = "TCGA-KIRC.GDC_phenotype_train.csv")

####
library(tidyverse)
surv = read.table(file = 'TCGA-KIRC.GDC_phenotype.tsv', sep = '\t', header = TRUE) 
surv <- TCGA.KIRC.GDC_phenotype
#整理生存信息数据
surv$submitter_id.samples <- substring(surv$submitter_id.samples,1,15) %>% gsub("-",".",.)
surv <- surv[!duplicated(surv$submitter_id.samples),]   #去重复

rownames(surv) <- surv$submitter_id.samples

comgene <- intersect(rownames(surv),rownames(ucox_surv.expr_01))#取交集

table(substr(comgene,14,15))#分组计数

surv <- surv[comgene,]

comgene <- intersect(rownames(surv),rownames(high))#取交集
surv <- surv[comgene,]
table(substr(comgene,14,15))#分组计数
write.csv(surv,file = "high.csv")

comgene <- intersect(rownames(surv),rownames(low))#取交集
surv <- surv[comgene,]
table(substr(comgene,14,15))#分组计数
write.csv(surv,file = "low.csv")




