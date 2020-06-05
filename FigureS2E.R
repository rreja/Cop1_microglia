## Heatmaps comparing select homeostatic and DAM genes
library(pheatmap)
eset.NGS1484 = readRDS('data/eset_NGS1484_rpkm.rds')
eset.NGS1484 = eset.NGS1484[,c("SAM24335239","SAM24335240","SAM24335241","SAM24335242","SAM24335243","SAM24335244","SAM24335245","SAM24335246")]

homeo = c('Tgfbr1','Mertk','Entpd1','P2ry12','Tmem119','Hexb','Smad3')
dams = c('Apoe','Itgax','Ctsb','Ctsd','Tyrobp','B2m','Fth1','Lyz2','Clec7a')

hRPKM = eset.NGS1484[which(fData(eset.NGS1484)$symbol %in% homeo),]
nRPKM = eset.NGS1484[which(fData(eset.NGS1484)$symbol %in% dams),]


dat = apply(exprs(hRPKM), 2, function(x) { x[which(is.na(x))] = 0; x } )
dat1 = apply(exprs(nRPKM), 2, function(x) { x[which(is.na(x))] = 0; x } )
dat = t(scale(t(dat)))
dat1 = t(scale(t(dat1)))
rownames(dat) = fData(hRPKM)$symbol
rownames(dat1) = fData(nRPKM)$symbol
dat = dat[apply(dat, 1, function(x) { !all(is.na(x))}),]
dat1 = dat1[apply(dat1, 1, function(x) { !all(is.na(x))}),]

## define col annotation
colAnnot = data.frame("Group"=hRPKM$genotype)
row.names(colAnnot) = colnames(dat)
colAnnot1 = data.frame("Group"=nRPKM$genotype)
row.names(colAnnot1) = colnames(dat1)


## Plot
pdf('figures/FigureS2E.pdf')
pheatmap(dat, annotation_col = colAnnot, cluster_cols=FALSE, border_color = NA, 
         cluster_rows = F, show_colnames = FALSE, fontsize_row = 10, show_rownames = T)

pheatmap(dat1, annotation_col = colAnnot1, cluster_cols=FALSE, border_color = NA, 
         cluster_rows = F, show_colnames = FALSE, fontsize_row = 10, show_rownames = T)

dev.off()