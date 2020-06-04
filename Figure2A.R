library(dplyr)
library(lattice)

## 4-way plot between LPS induced change vs changes due to Cop1 deletion
lfc_LPS = readRDS('data/NGS1307_LPS_treatment.rds')
lfc_cop1d = readRDS('data/NGS1307_Cop1_deletion.rds')
colnames(lfc_cop1d) = c("ID","symbol","lfc.wt.cop1d","p.wt.cop1d")

## Merge the two datasets
merged = left_join(lfc_LPS,lfc_cop1d,by='ID')
## Create 4 - way plot
red = "#B3252E"; redfill = "#e6550d"
green = "#1D9062"; greenfill = "#31a354"
blue = "#1361AD"; bluefill = "#068ae5"
grey = "#d6d4d4" ; greyfill="#d6d4d4"
black = "#000000"; blackfill = "#000000"

col = rep(grey, nrow(merged))
# Both up
col[which(merged$lfc.wt.lps.untreated > 1 & merged$p.wt.lps.untreated < 0.05 & merged$lfc.wt.cop1d > 1 & merged$p.wt.cop1d < 0.05)] = red
# Both down
col[which(merged$lfc.wt.lps.untreated < -1 & merged$p.wt.lps.untreated < 0.05 & merged$lfc.wt.cop1d < -1 & merged$p.wt.cop1d < 0.05)] = blue
# Up in only one or the other
col[which(merged$lfc.wt.lps.untreated > 1 & merged$p.wt.lps.untreated < 0.05 & merged$lfc.wt.cop1d < -1 & merged$p.wt.cop1d < 0.05)] = green
col[which(merged$lfc.wt.lps.untreated < -1 & merged$p.wt.lps.untreated < 0.05 & merged$lfc.wt.cop1d > 1 & merged$p.wt.cop1d < 0.05)] = green


## To add regression line
y = lm(merged$lfc.wt.cop1d  ~ merged$lfc.wt.lps.untreated, data = merged)
## Add corr coeff
corr = cor.test(merged$lfc.wt.lps.untreated,merged$lfc.wt.cop1d)

# Plotting

pch = rep(19, nrow(merged))
cex = rep(0.6, nrow(merged))

cex[which(merged$lfc.wt.lps.untreated > 1 & merged$p.wt.lps.untreated < 0.05 & merged$lfc.wt.cop1d > 1 & merged$p.wt.cop1d < 0.05)] = 0.8

## Label following genes
glist = c('Saa3','Fpr1','Fpr2','Mmp3','Mx1','Cxcl2','Ifi205','Cxcl10','Ccl5','Irg1','Cd69','Slfn4')
df = merged[merged$symbol.x %in% glist,]

pdf('figures/Figure2A.pdf')
xyplot(merged$lfc.wt.cop1d ~ merged$lfc.wt.lps.untreated, pch=pch, ylab=" Changes induced by COP1 deletion",
       xlab="LPS induced changes in COP1 wild-type microglia", main = paste0('R = ',round(corr$estimate,2)),
       par.settings = list(plot.symbol = list(col = col,cex = cex,pch=pch)),
       panel=function(...) { panel.xyplot(...)
         panel.abline(h=1, lty = "dashed", col = "#63636388")
         panel.abline(h=-1, lty = "dashed", col = "#63636388")
         panel.abline(v=1, lty = "dashed", col = "#63636388")
         panel.abline(v=-1, lty = "dashed", col = "#63636388")
         panel.abline(h=0, lty = "solid", col = "#63636388")
         panel.abline(v=0, lty = "solid", col = "#63636388")
         panel.abline(a = 0, b = 1, lty = "solid", col = "#63636388")
         panel.text(x= df$lfc.wt.lps.untreated, y = df$lfc.wt.cop1d, labels = df$symbol.x, col = 'black', pos=3)
       })

dev.off()