####差异分析####
KIRC_counts_all <- TCGA_KIRC_counts_log2+1
comgene <- intersect(rownames(Diffsing.Corr.output),rownames(KIRC_counts_all))#取交集
table(substr(comgene,14,16))
pearson_gene_exp <- KIRC_counts_all[comgene,]
write.csv(pearson_gene_exp,file = "pearson_gene_exp_counts.csv")

library(tidyverse)
#安装BiocManager
if(!require(DESeq2))BiocManager::install('DESeq2')#DESeq2是用来跑差异分析的
library(DESeq2)
#counts为Pearson找出的差异基因的表达谱
counts <- pearson_gene_exp
counts <- ceiling(2^(counts)-1)

tumor <- colnames(counts)[substr(colnames(counts),14,15) == "01"]
counts_01A <- counts[,tumor]
normal <- colnames(counts)[substr(colnames(counts),14,15) == "11"]
counts_11A <- counts[,normal]
counts <- cbind(counts_11A,counts_01A)

counts = counts[apply(counts, 1, function(x) sum(x > 1) > 32), ]
conditions=data.frame(sample=colnames(counts),
                      group=factor(ifelse(substr(colnames(counts),14,15) == "01","T","N"),levels = c("N","T"))) %>% 
  column_to_rownames("sample")

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = conditions,
  design = ~ group)
dds <- DESeq(dds)#跑差异分析的代码
resultsNames(dds)
res <- results(dds)
save(res,file = "KIRC_DEG_result.rda")#一定要保存！
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 3, padj < 0.05)#根据自己需要#表示调节P值小于0.05，癌症基因表达相对正常翻倍的倍数的绝对值增大三倍
write.csv(res_deseq2, file = "res_deseq2.csv")
#DEG:differentially expressed genes
#DEG(LUSC_DEG_result)是计算出来的差异分析的结果，res_deseq2是根据一定标准的筛选结果
#读取差异基因文件
#直接在文件夹双击

####TCGA差异分析热图####
setwd("xena")
library(tidyverse)
#exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#DEG就是差异基因的表达结果，在文件里是LUSC_DEG_result.rda（res）
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 3, padj < 0.05)

logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)


library(pheatmap)
cg = rownames(DEG)[DEG$change !="NOT"]
exp <- counts#从上面的代码可以知道，这里的exp也是差异分析基因的表达谱
exp_diff <- exp[cg,]
group_list=factor(ifelse(substr(colnames(exp),14,15) == "01","T","N"),levels = c("N","T"))

annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(exp_diff)
pheatmap(exp_diff,
         annotation_col=annotation_col,#需要的文件有表达谱exp-diff,分组信息group-list
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
dev.off()

####TCGA差异分析火山图####
#差异分析火山图和热图的数据是一样的
setwd("xena")
library(tidyverse)
#exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 0, padj < 0.05)

logFC_cutoff <- 3
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)


#install.packages("ggpubr")
#install.packages("ggthemes")
library(ggpubr)
library(ggthemes)

DEG$logP <- -log10(DEG$padj)
ggscatter(DEG,
          x = "log2FoldChange", y = "logP") +
  theme_base()#开始画图

#增加基因上下调信息
ggscatter(DEG, x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1) +
  theme_base()

#添加分界线
ggscatter(DEG, x = "log2FoldChange", y = "logP", xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1) +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-3, 3), linetype = "dashed")#其中的（-1，1)这个范围可改变
#dev.off()

#添加基因标签信息
DEG$Label = ""   #新加一列label
DEG <- DEG[order(DEG$padj), ]   #对差异基因的p值进行从小到大的排序
DEG$Gene <- rownames(DEG)#DEG里新增一列叫做Gene
#高表达的基因中，选择fdr值最小的5个
up.genes <- head(DEG$Gene[which(DEG$change == "UP")], 5)
#低表达的基因中，选择fdr值最小的5个
down.genes <- head(DEG$Gene[which(DEG$change == "DOWN")], 5)
#将up.genes和down.genes合并，并加入到Label中
DEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
DEG$Label[match(DEG.top5.genes, DEG$Gene)] <- DEG.top5.genes

ggscatter(DEG, x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1,
          label = DEG$Label,
          font.label = 8,
          repel = T,
          xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)") +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-3, 3), linetype = "dashed")

dev.off()
