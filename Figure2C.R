# Let's make heatmap visualizaion for Brad's modules
library(pheatmap)
library(dplyr)
eset.NGS1484 = readRDS('data/eset_NGS1484_rpkm.rds')
## Arrange the COP1 WT and KO replicates accordingly
eset.NGS1484 = eset.NGS1484[,c("SAM24335239","SAM24335240","SAM24335241","SAM24335242","SAM24335243","SAM24335244","SAM24335245","SAM24335246")]

tt = read.table('data/Brads_modules_with_gene_symbol.txt',sep="\t",header = T, stringsAsFactors = F)
mRPKM = eset.NGS1484[which(fData(eset.NGS1484)$symbol %in% tt$symbol),]

dat = apply(exprs(mRPKM), 2, function(x) { x[which(is.na(x))] = 0; x } )
dat = t(scale(t(dat)))
rownames(dat) = fData(mRPKM)$symbol
dat = dat[match(tt$symbol, rownames(dat)),]
dat = dat[apply(dat, 1, function(x) { !all(is.na(x))}),]

## define col annotation
colAnnot = data.frame("Group"=mRPKM$genotype)
row.names(colAnnot) = colnames(dat)

## define row annotation
annotation_row = data.frame( Geneset = tt$geneSet)
row.names(annotation_row) = rownames(dat)

pdf('figures/Figure2C.pdf')
## Plot
pheatmap(dat, annotation_col = colAnnot, annotation_row = annotation_row, cluster_cols=FALSE, border_color = NA, 
             cluster_rows = F, show_colnames = FALSE, fontsize_row = 6,treeheight_row = 0, show_rownames = T)
dev.off()
