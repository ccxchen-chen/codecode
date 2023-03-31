library(GSVA)
library(limma)
library(GSEABase)

#inputFile="symbol.txt"    
gmtFile="immune.gmt"    
      

rt=`TCGA_KIRC_fpkm_log2+1`
#rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
#rownames(rt)=rt[,1]
#exp=rt[,2:ncol(rt)]
#dimnames=list(rownames(exp),colnames(exp))
mat <- rt#mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]
geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}

ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut,file="ssgseaOut.txt",sep="\t",quote=F,col.names=F)


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056
