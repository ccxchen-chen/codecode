#多因素coxgene可视化热图
#读取risk_train,在读进来之前用excel给高低风险排序
counts <- risk_all
counts <- counts[,-1]
counts <- counts[,-1]
counts <- counts[,-1]
counts1 <- counts[,-ncol(counts)] 
counts1 <- counts1[,-ncol(counts1)]
counts1 <- t(counts1)


library(pheatmap) # 加载pheatmap这个R包
#group_list=factor(ifelse(substr(counts$risk) == "high","high","low"),levels = c("low","high"))
group_list=counts$risk
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=rownames(counts)
#mydata1[,c(3:7)] <- as.numeric(unlist(mydata1[,c(3:7)]))

pheatmap(counts1,
         annotation_col=annotation_col,#需要的文件有表达谱exp-diff,分组信息group-list
         annotation_colors = list(group = c("low" ="#01468b","high"= "#ee0000")),
         scale = "row",
         show_rownames = T,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         fontsize = 10,
         fontsize_row=8,
         fontsize_col=3)

