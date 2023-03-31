dat <- phenotype_sur_score


library(ggpubr)  
p <- ggboxplot(dat, x = "race", y = "riskScore",#x代表生存状态
               color = "race", palette = "jco",
               add = "jitter")
p <- p + stat_compare_means()   #  Add p-value
p   #得出预测结果
dev.off()
