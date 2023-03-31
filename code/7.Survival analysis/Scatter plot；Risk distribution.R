
###5. 绘制散点分布图
#install.packages("ggpubr")

dat <- risk_all


library(ggpubr)  
p <- ggboxplot(dat, x = "fustat", y = "riskScore",#x代表生存状态
               color = "fustat", palette = "jco",
               add = "jitter")
p <- p + stat_compare_means()   #  Add p-value
p   #得出预测结果


#画风险分布图
#y=生存时间
rt <- dat#这里进行赋予只是因为下面的代码都是rt
color=as.vector(rt$fustat)
color[color==1]="indianred1"
color[color==0]="lightseagreen"
plot(rt$futime, pch=19,
     xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("indianred1","lightseagreen"),cex=1.2)
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
abline(v=lowLength,lty=2)
dev.off()

#y=riskscore
rt <- rt[order(rt[,12]),]#通过rt数据框的第25列（即风险评分）进行排序（此处默认升序，如果想降序排列，在逗号后面加decreasing=t）
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
plot(line, type="p", pch=20,
     xlab="Patients (increasing risk socre)", ylab="Risk score",
     col=c(rep("lightseagreen",lowLength),rep("indianred1",highLength)) )
abline(h=median(rt$riskScore),v=lowLength,lty=2)#风险评分的中位数、低组长度
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("indianred1","lightseagreen"),cex=1.2)
dev.off()
