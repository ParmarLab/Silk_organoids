################################################################################
# Silk scaffolding drives self-assembly of functional and mature human brain organoids
################################################################################

################################################################################
# Data preparation 
################################################################################

# Packages
library(DESeq2)
library(ggplot2)
library(cowplot)
library(EnhancedVolcano)
library("AnnotationDbi")
library("org.Hs.eg.db")
library("pheatmap")
library("RColorBrewer")

load("../Silk scaffolding drives self-assembly of functional and mature human brain organoids - RNAseq Data/Bulk RNA seq/rnaseq_fb_org.Rdata")
res <- results(dds)
res <-res[complete.cases(res),]
res.s <- res[res$padj<0.05,]
res <- lfcShrink(dds, coef="treat_std_vs_silk", type="apeglm")
ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

res$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
resOrdered <- res[order(res$padj),]
write.csv(resOrdered,file = "results.deseq.csv")

################################################################################
# Figure 2 
################################################################################

# Figure 2 A
################################################################################        
hm.genes <- read.csv("../Silk scaffolding drives self-assembly of functional and mature human brain organoids - RNAseq Data/Bulk RNA seq/Gene_list_heatmap_Fig2A.csv",header = T, sep = ";",
                     skip = 1)
colnames(hm.genes)<-c("Category","Ensembl","Symbol")
hm.e <- assay(dds)[hm.genes$Ensembl,]

rownames(hm.e)<-hm.genes$Symbol
rowAnno <- data.frame(Category=hm.genes[,c("Category")],row.names = hm.genes$Symbol)

figure2A <- pheatmap(hm.e,scale = "row",
            annotation_row = rowAnno[],
            cluster_rows = F,border_color = "darkgrey",
            color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100))
figure2A

################################################################################
# Supplementary Figure 1 
################################################################################

# Supplementary Figure 1 G
################################################################################
keyvals <- ifelse(
  res$log2FoldChange < -1, "#EECF82",
  ifelse(res$log2FoldChange > 1, "#7AAABA",
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == "#7AAABA"] <- 'Conventional'
names(keyvals)[keyvals == 'black'] <- 'Not DEG'
names(keyvals)[keyvals == "#EECF82"] <- 'Silk'

suppfigure1G <- EnhancedVolcano(res,
                lab=res$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c("RSPO2","RSPO3","RSPO1","LMX1A","HEY2","BEGAIN",
                              "CCN2","DMRT3","AFP","HOXB1","GATA6","WNT1","TBX3",
                              "DLK1","FZD5","SIX3","FGF8","EPHA2","PAX8","ATP1A2",
                              "SOX17","OPTC"),
                xlab = bquote(~Log[2]~ 'fold change'),
                labSize = 7,
                colCustom = keyvals,
                colAlpha = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                typeConnectors = "closed",
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1,
                borderColour = 'black')
suppfigure1G

sessionInfo()
#R version 4.1.0 (2021-05-18)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Big Sur 10.16

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#  [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] RColorBrewer_1.1-2          pheatmap_1.0.12             org.Hs.eg.db_3.13.0        
#[4] AnnotationDbi_1.54.1        EnhancedVolcano_1.10.0      ggrepel_0.9.1              
#[7] cowplot_1.1.1               ggplot2_3.3.5               DESeq2_1.32.0              
#[10] SummarizedExperiment_1.22.0 Biobase_2.52.0              MatrixGenerics_1.4.3       
#[13] matrixStats_0.61.0          GenomicRanges_1.44.0        GenomeInfoDb_1.28.4        
#[16] IRanges_2.26.0              S4Vectors_0.30.2            BiocGenerics_0.38.0  

rm(list = ls())