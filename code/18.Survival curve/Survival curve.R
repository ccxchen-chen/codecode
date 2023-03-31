

comgene <- intersect(rownames(risk_all),rownames(phenotype_sur_score))#取交集
table(substr(comgene,14,16))
risk_all <- risk_all[comgene,]
phenotype_sur_score <- phenotype_sur_score[comgene,]
data <- cbind(phenotype_sur_score,risk_all$risk)

colnames(data)[6] <- 'risk'
#通过这句代码筛选不同条件的临床信息,不能只有等于号
data1 <- subset(data,age>65)
####生存分析####
surv <- data1

surv$futime <- surv$futime/12#把生存时间改为以月为单位

surv$risk <- factor(surv$risk, levels = c("low","high")) #把low、high改为一字数字类型
class(surv$risk)
table(surv$risk)

#计算HR
library(survival)
data.survdiff <- survdiff(Surv(futime, fustat) ~ risk,data      = surv)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))



#install.packages("survival")
library(survival)
fitd <- survdiff(Surv(futime, fustat) ~ risk,#进行生存分析
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线
fit <- survfit(Surv(futime, fustat)~ risk, data = surv)
summary(fit)
###添加P值
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
text(25, 0.2, p.lab)
#3. 绘制生存曲线

#方法2
#install.packages("survminer")
library(survminer)
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表（可以去掉）
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,120), # x轴长度，一般为0-10年
           break.time.by = 20, # x轴步长为20个月
           legend.title = "",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Months)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块（可以去掉）
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()

