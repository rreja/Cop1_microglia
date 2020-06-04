library(Seurat)
library(ggplot2)
## Single cell RNAseq figures
obj = readRDS('data/Seurat_obj.rds')

pdf('figures/Figure2D.pdf')
# cluster plot
DimPlot(object = obj, cols = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#cab2d6'))

## genotype plot
DimPlot(obj, group.by = 'genotype', cols = c('#ee3424','#000000'))

## Plotting markers
p1 <- FeaturePlot(obj,c('Apoe','P2ry12'), combine = FALSE)
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'blue'),  limits = c(0, 8))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)

dev.off()