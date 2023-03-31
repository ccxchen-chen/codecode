
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")
TCGA_CIBERSORT_Results_fromRcode<- TCGA_CIBERSORT_Results_fromRcode[!duplicated(TCGA_CIBERSORT_Results_fromRcode$X),]   #去重复
rownames(TCGA_CIBERSORT_Results_fromRcode) <- TCGA_CIBERSORT_Results_fromRcode[,1]
com <- intersect(rownames(risk_all),rownames(TCGA_CIBERSORT_Results_fromRcode))#列名取交集，为了使两个表的列名变得相同
risk_all <- risk_all[com,]
TCGA_CIBERSORT_Results_fromRcode <- TCGA_CIBERSORT_Results_fromRcode[com,]
TCGA_CIBERSORT_Results_fromRcode$riskScore <- risk_all$riskScore
#加载包
library(ggplot2)
library(ggpubr)
library(ggExtra)

inputFile="input.txt"      
gene1="riskScore"             #第一个基因名字
gene2="Neutrophils"              #第二个基因名字
   

#读取输入文件，提取基因表达量
rt <- TCGA_CIBERSORT_Results_fromRcode
#rt=read.table(inputFile,sep="\t",header=T,check.names=F,row.names=1)
x=as.numeric(rt[,gene1])
y=as.numeric(rt[,gene2])

#相关性分析
df1=as.data.frame(cbind(x,y))
corT=cor.test(x,y,method="spearman")
cor=corT$estimate
pValue=corT$p.value
p1=ggplot(df1, aes(x, y)) + 
			xlab(gene1)+ylab(gene2)+
			geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
			stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))

#出图
pdf(file="Neutrophils cor.pdf",width=5,height=4.8)
print(p1)
dev.off()

#出图2
pdf(file="Neutrophils cor.density.pdf",width=5,height=5)
print(p2)
dev.off()

