rm(list = ls())
require(maftools)
options(stringsAsFactors = F) 
library(data.table)
tmp <- fread("TCGA-KIRC.mutect2_snv.tsv.gz")
head(tmp)   
colnames(tmp) =c( "Tumor_Sample_Barcode", "Hugo_Symbol", 
                  "Chromosome", "Start_Position", 
                  "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                  "HGVSp_Short" , 'effect' ,"Consequence",
                  "vaf" ) 
tmp$Entrez_Gene_Id =1
tmp$Center ='ucsc'
tmp$NCBI_Build ='GRCh38'
tmp$NCBI_Build ='GRCh38'
tmp$Strand ='+'
tmp$Variant_Classification = tmp$effect
tail(sort(table(tmp$Variant_Classification )))
tmp$Tumor_Seq_Allele1 = tmp$Reference_Allele
tmp$Variant_Type = ifelse(
  tmp$Reference_Allele %in% c('A','C','T','G') & tmp$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
  'SNP','INDEL'
)
table(tmp$Variant_Type )

tmp$Tumor_Sample_Barcode <- substr(tmp$Tumor_Sample_Barcode,1,12)









tcga.coad = read.maf(maf = tmp,
                     vc_nonSyn=names(tail(sort(table(tmp$Variant_Classification )))))

plotmafSummary(maf = tcga.coad, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
oncoplot(maf = tcga.coad, top = 10) # 高频突变的前10个基因
