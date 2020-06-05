# heatmap visualizaion for Brad's modules for aged mice
library(pheatmap)
library(dplyr)
eset.NGS2450 = readRDS('data/eset.NGS2450.rpkm.rds')
## Arrange the COP1 WT and KO replicates accordingly
eset.NGS2450 = eset.NGS2450[,c("SAM24355936","SAM24355937","SAM24355938","SAM24355939","SAM24355940","SAM24355941")]

tt = read.table('data/Brads_modules_with_gene_symbol.txt',sep="\t",header = T, stringsAsFactors = F)
mRPKM = eset.NGS2450[which(fData(eset.NGS2450)$symbol %in% tt$symbol),]

dat = apply(exprs(mRPKM), 2, function(x) { x[which(is.na(x))] = 0; x } )
dat = t(scale(t(dat)))
rownames(dat) = fData(mRPKM)$symbol
dat = dat[match(tt$symbol, rownames(dat)),]
dat = dat[apply(dat, 1, function(x) { !all(is.na(x))}),]

## define col annotation
colAnnot = data.frame("Group"=mRPKM$group)
row.names(colAnnot) = colnames(dat)

## define row annotation
tt = tt[tt$symbol %in% rownames(dat),]
annotation_row = data.frame( Geneset = tt$geneSet)
row.names(annotation_row) = rownames(dat)

pdf('figures/FigureS2D.pdf')
## Plot
pheatmap(dat, annotation_col = colAnnot, annotation_row = annotation_row, cluster_cols=FALSE, border_color = NA, 
         cluster_rows = F, show_colnames = FALSE, fontsize_row = 6,treeheight_row = 0, show_rownames = T)
dev.off()
