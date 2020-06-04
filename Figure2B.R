## Volcano plot for Cop1 deletion in vivo
library(EnhancedVolcano)
lfc_1484 = readRDS('data/NGS1484_Cop1_deletion.rds')

# Adjust point size
lfc_1484$ptsize = 0.8
lfc_1484$ptsize <- ifelse(((lfc_1484$lfc.1484 >1 |lfc_1484$lfc.1484 < -1) & (lfc_1484$p.1484 < 0.05)) , 2, 0.8)

pdf('figures/Figure2B.pdf')
EnhancedVolcano(lfc_1484,lab = lfc_1484$symbol, x = 'lfc.1484', y = 'p.1484', 
                xlim = c(-5.5, 10), selectLab = c('Apoe','C3','Fth1','Cxcr2','Saa3','Il12b',
                                                  'P2ry12','Cx3cr1','Itgax','Il12b','Ifitm1',
                                                  'Ifi204','Cxcl10','Ccl5','Cxcl2','Gpnmb','Irg1'),
                DrawConnectors = T, 
                widthConnectors = 1.0, 
                colConnectors = 'black',
                gridlines.major = F,
                gridlines.minor = F,
                FCcutoff = 1,
                col = c("grey30", "grey30", "royalblue", "red2"),
                transcriptPointSize = lfc_1484$ptsize,
                colAlpha = 1.0)

dev.off()