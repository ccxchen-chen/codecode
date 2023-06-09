

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("scales")
#install.packages("ggplot2")
#install.packages("ggtext")
#install.packages("tidyverse")
#install.packages("ggpubr")


#???ð?
library(limma)
library(scales)
library(ggplot2)
library(ggtext)
library(reshape2)
library(tidyverse)
library(ggpubr)

riskFile="risk_all.txt"      #?????ļ?
immFile="timer_risk.csv"     #????ϸ???????ļ?
#setwd("C:\\biowolf\\Necroptosis\\25.immuneCor")     #???ù???Ŀ¼

#??ȡ?????????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#??ȡ????ϸ???????ļ?
immune=read.csv(immFile, header=T, sep=",", check.names=F, row.names=1)
immune=as.matrix(immune)
rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))
immune=avereps(immune)

#对风险文件和免疫细胞浸润文件取交集，得到交集样品
sameSample=intersect(row.names(risk), row.names(immune))
risk=risk[sameSample, "riskScore"]
immune=immune[sameSample,]

#对风险得分和免疫细胞进行相关性分析
x=as.numeric(risk)
x[x>quantile(x,0.99)]=quantile(x,0.99)
outTab=data.frame()
for(i in colnames(immune)){
	y=as.numeric(immune[,i])
	if(sd(y)<0.001){next}
	corT=cor.test(x, y, method="spearman")
	cor=corT$estimate
	pvalue=corT$p.value
	if(pvalue<0.05){
		outTab=rbind(outTab,cbind(immune=i, cor, pvalue))
		#??????????ɢ??ͼ
		outFile=paste0("cor.", i, ".pdf")
		outFile=gsub("/", "_", outFile)
		df1=as.data.frame(cbind(x,y))
		p1=ggplot(df1, aes(x, y)) + 
				  xlab("Risk score") + ylab(i)+
				  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
				  stat_cor(method = 'spearman', aes(x =x, y =y))
		#??????ͼ??
		pdf(file=outFile, width=5, height=4.7)
		print(p1)
		dev.off()
	}
}
#?????????Խ???
write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)

#????????ͼ
corResult=read.table("corResult.txt", head=T, sep="\t")
corResult$Software=sapply(strsplit(corResult[,1],"_"), '[', 2)
corResult$Software=factor(corResult$Software,level=as.character(unique(corResult$Software[rev(order(as.character(corResult$Software)))])))
b=corResult[order(corResult$Software),]
b$immune=factor(b$immune,levels=rev(as.character(b$immune)))
colslabels=rep(hue_pal()(length(levels(b$Software))),table(b$Software))     #??????ɫ
#????ͼ??
pdf(file="correlation.pdf", width=9, height=8)
ggplot(data=b, aes(x=cor, y=immune, color=Software))+
	labs(x="Correlation coefficient",y="Immune cell")+
	geom_point(size=4.1)+
	theme(panel.background=element_rect(fill="white",size=1,color="black"),
	      panel.grid=element_line(color="grey75",size=0.5),
	      axis.ticks = element_line(size=0.5),
	      axis.text.y = ggtext::element_markdown(colour=rev(colslabels)))
dev.off()
