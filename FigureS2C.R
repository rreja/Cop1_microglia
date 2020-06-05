## 4-way plot comparing the COP1 KO changes in young and old mice
library(lattice)

merged = read.table('data/Merged_lfcPval_young_old_mice.txt',header = T,stringsAsFactors = F,sep = "\t")
## Plotting colors order: grey -> green -> red -> blue
merged$order <- ifelse(((merged$lfc.1484 > 1 | merged$lfc.1484 < -1) & (merged$p.1484 < 0.05 & merged$p.2450 > 0.05)),2,1)
merged[((merged$lfc.2450 > 1 | merged$lfc.2450 < -1) & (merged$p.2450 < 0.05 & merged$p.1484 > 0.05)), 'order'] <- 3
merged[(merged$lfc.1484 > 1 & merged$p.1484 < 0.05 & merged$lfc.2450 > 1 & merged$p.2450 < 0.05), 'order'] <- 4
merged[(merged$lfc.1484 < -1 & merged$p.1484 < 0.05 & merged$lfc.2450 < -1 & merged$p.2450 < 0.05), 'order'] <- 4
merged = merged[order(merged$order),]


## Create 4 - way plot
red = "#e6550d"; 
green = "#31a354"; 
blue = "#068ae5"; 
grey = "#d6d4d4" ; 


col = rep(grey, nrow(merged))
# Both up
col[which(merged$lfc.1484 > 1 & merged$p.1484 < 0.05 & merged$lfc.2450 > 1 & merged$p.2450 < 0.05)] = blue
# Both down
col[which(merged$lfc.1484 < -1 & merged$p.1484 < 0.05 & merged$lfc.2450 < -1 & merged$p.2450 < 0.05)] = blue
# Up in only one or the other
col[which((merged$lfc.1484 > 1 | merged$lfc.1484 < -1) & (merged$p.1484 < 0.05 & merged$p.2450 > 0.05))] = green
col[which((merged$lfc.2450 > 1 | merged$lfc.2450 < -1) & (merged$p.2450 < 0.05 & merged$p.1484 > 0.05))] = red


## To add regression line
y = lm(merged$lfc.1484  ~ merged$lfc.2450, data = merged)
## Add corr coeff
corr = cor.test(merged$lfc.2450,merged$lfc.1484)

# Plotting

pch = rep(19, nrow(merged))
cex = rep(0.4, nrow(merged))

cex[which(merged$lfc.1484 > 1 & merged$p.1484 < 0.05 & merged$lfc.2450 > 1 & merged$p.2450 < 0.05)] = 1
cex[which(merged$lfc.1484 < -1 & merged$p.1484 < 0.05 & merged$lfc.2450 < -1 & merged$p.2450 < 0.05)] = 1

## Label following genes
glist = c('Cxcr2','Ifitm1','Ctse','Apoe','Cxcl2','Ifi205',
          'Cxcl10','Itgax','Ccl5','C3','Il12b','Gpnmb',
          'Irg1','P2ry12','P2ry13','Cx3cr1','Il12a','Hist1h1d',
          'Bmx','Tox2','Hist1h2ad','Mfap4','Rnasel','Zdbf2','Evc2',
          'Ifnb1','Dupd1','Saa3','Il12b','Ifi205','Cd69','Jaml')
df = merged[merged$symbol.x %in% glist,]

pdf('figures/FigureS2C.pdf')
xyplot(merged$lfc.1484 ~ merged$lfc.2450, pch=pch, ylab="Log FC (2-3 weeks post Cop1 deletion)",
       xlab="Log FC (>1.5 years post Cop1 deletion)", main = paste0('R = ',round(corr$estimate,2)),
       par.settings = list(plot.symbol = list(col = col,cex = cex,pch=pch)),
       panel=function(...) { panel.xyplot(...)
         panel.abline(h=1, lty = "dashed", col = "#63636388")
         panel.abline(h=-1, lty = "dashed", col = "#63636388")
         panel.abline(v=1, lty = "dashed", col = "#63636388")
         panel.abline(v=-1, lty = "dashed", col = "#63636388")
         panel.abline(h=0, lty = "solid", col = "#63636388")
         panel.abline(v=0, lty = "solid", col = "#63636388")
         panel.abline(a = 0, b = 1, lty = "solid", col = "#63636388")
         panel.text(x= df$lfc.2450, y = df$lfc.1484, labels = df$symbol.x, col = 'black', pos=3)
       })

dev.off()