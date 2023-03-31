#install.packages("ggpubr")
#去TIDE官网输入全部fpkm.txt文件

library(ggpubr)     
Group=risk_all       
score=TIDE_result      

#数据合并
sameSample=intersect(row.names(Group), row.names(score))
data=cbind(score[sameSample,,drop=F], Group[sameSample,"risk",drop=F])
data[,"Group"]=factor(data[,"risk"], levels=c("high","low"))

#设置比较组
group=levels(factor(data$risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制小提琴图
bioCol=c("#FF0000","#0066FF","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
for(i in colnames(data)[1:(ncol(data)-1)]){
	violin=ggviolin(data, x="risk", y=i, fill = "risk",
	         xlab="risk", ylab=i,
	         legend.title="risk",
	         palette=bioCol, 
	         add="boxplot", add.params = list(fill="white"))+ 
		stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	pdf(file=paste0("vioplot.", i, ".pdf"), width=6, height=5.5)
	print(violin)
	dev.off()
}


