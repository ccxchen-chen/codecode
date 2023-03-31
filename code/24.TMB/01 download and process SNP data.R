###################################################################

############安装TCGAbiolinks包#######################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
BiocManager::install("maftools")

##############加载TCGAbiolinks包                  
library(TCGAbiolinks)

##############利用GDCquery函数下载SNP原始数据
query_SNV <- GDCquery(project = "TCGA-KIRC",
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(query_SNV)
save(query_SNV,file = "data/01 query_SNV.Rdata")

##############提取原始数据得到矩阵
snv <- GDCprepare(query_SNV)
save(snv,file = "snv.Rdata")

###############读取临床分组文件
library(data.table)
clin<- risk_all
clincluster <- clin[,c("id","risk")]
clincluster$id <- gsub("\\.","-",clincluster$id)

snv$Tumor_Sample_Barcode <- substr(snv$Tumor_Sample_Barcode,1,15)

snv.high <- snv[(snv$Tumor_Sample_Barcode %in% clincluster$id[clincluster$risk=="high"]),]
snv.low <- snv[(snv$Tumor_Sample_Barcode %in% clincluster$id[clincluster$risk=="low"]),]
save(snv.high,file = "snv.high.Rdata")
save(snv.high,file = "snv.low.Rdata")

#################读取数据为maf文件
library(maftools)
maf <- read.maf(maf = snv)
maf.high <- read.maf(maf = snv.high)
maf.low <- read.maf(maf = snv.low)

################瀑布图
oncoplot(maf = maf.high,
        top = 20)
oncoplot(maf = maf.low,
         top = 20)

################计算TMB
x = tmb(maf = maf,
          captureSize = 20,
          logScale = TRUE)   

head(x)
TMB <- x[,c(1,3)]


library(dplyr)
TMB_group <- mutate(TMB,"group")

for ( i in 1:575) {print(i)
  if(TMB_group[i,2]>2.25){TMB_group[i,3]="TMB-H"}
  if(TMB_group[i,2]<=2.25){TMB_group[i,3]="TMB-L"}
}
colnames(TMB_group) <- c("id","Tumor_tmbation_burden","group")
TMB_group <- as.data.frame(TMB_group)


colnames(TMB) <- c("id","Tumor_tmbation_burden")
TMB_cluster <- merge(TMB,clincluster,by="id")
TMB_cluster <- as.data.frame(TMB_cluster)
TMB_cluster$risk=ifelse(TMB_cluster$risk=="high", "High-risk", "Low-risk")
TMB_cluster$risk <- as.factor(TMB_cluster$risk)


TMB_compare <- merge(TMB_group,TMB_cluster,by="id")###重要
TMB_compare <- TMB_compare [,-4]
TMB_compare$risk <- as.factor(TMB_compare $risk)
TMB_compare$group <- as.factor(TMB_compare $group)
write.csv(TMB_compare,file = "TMB_compare.csv")
write.csv(TMB_group,file = "TMB_group.csv")


####################定义组间关系
group=levels(factor(TMB_cluster$risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

################ggplot2作图
library(ggplot2)
library(ggpubr)

ggplot(TMB_cluster, aes(x = risk, y = Tumor_tmbation_burden))+ #注释：”x=”，”y=”表示x轴和y轴表示的变量数值，p表示图像对象
geom_violin(aes(fill=risk),cex=1.2)+ #注释：画出violin plot的函数
  scale_fill_manual(values = c('#579ABB','#B978AE'))+ #定义颜色
  geom_boxplot(width=0.1,cex=1.2)+ #箱图及大小调整
  theme_classic(base_size = 15)+ 
  theme(axis.text = element_text(color = 'black'), #默认坐标轴刻度颜色其实是灰色，这里我们修改成黑色
        legend.position = 'none')+       #去掉图例

   stat_compare_means(comparisons = my_comparisons)

ggplot(TMB_compare,aes(x=risk,y=Tumor_tmbation_burden.x,fill=group))+
  geom_bar(stat="identity",position="fill")#+
  geom_text(aes(label=perc),position=position_stack(vjust=0.5))
  
save.image(file = "TMB.Rdata")
