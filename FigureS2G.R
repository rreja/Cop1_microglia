library(Seurat)
library(RColorBrewer)
## Heatmap for markers of cluster
obj = readRDS('data/Seurat_obj.rds')
all.genes = rownames(obj)
scaled_obj <- ScaleData(object = obj, verbose = FALSE,block.size = 2000, features = all.genes)
genes = c('P2ry12','P2ry13','Tmem119','Cx3cr1','Cst3','Gpr34','Sall1','Siglech','Sparc','Arhgap5',
          'Apoe','Capg','Milr1','Ctsb','Ctsz','Cd34',
          'Irf7','Isg15','Ifi27l2a','Ifi204','Sp100','Ifit3','Rtp4','Xaf1',
          'Ms4a7','Mrc1','Pf4','Dab2','F13a1','F13a1',
          'Ube2c','Birc5','Top2a','Ccnb2','Ccnb2','Nusap1','Cenpf','Cdca8','Mki67','Cdca3',
          'Trem3','Clec4e','Hp','Ifitm6','Gda','S100a6','Plbd1','Msrb1','Flna','Itgal')

breaklist = seq(-2, 2, by = 0.001)
color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaklist))
DoHeatmap(scaled_obj, features = genes) + scale_fill_gradientn(colors = color)
