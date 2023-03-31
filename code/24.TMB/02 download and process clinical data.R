############################################################
#################清理环境变量##################
rm(list = ls())

################准备工作####################
#######从TCGA中下载临床数据################

########处理临床数据#################
clin <- read.table("resource/clinical.cart.2022-11-11/clinical.tsv",header = T, quote="", sep = "\t")
clin1 <- clin[,c("case_id","case_submitter_id","vital_status","days_to_last_follow_up","ajcc_pathologic_t","ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stage")]
colnames(clin1) <- c("case_id","case_submitter_id","vital_status","days_to_last_follow_up","ajcc_pathologic_t","ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stage")
clin1 <- clin1[clin1$days_to_last_follow_up!="--",]
clin1$days_to_last_follow_up <- as.numeric(clin1$days_to_last_follow_up)
clin1 <- clin1[clin1$vital_status!="--",]
clin1$vital_status <- ifelse(clin1$vital_status=="Alive",0,1)
