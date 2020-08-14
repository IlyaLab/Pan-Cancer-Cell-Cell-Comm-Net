#cell 2 cell analysis

library(dplyr)
library(tidyr)
library(ggridges)
library(stringr)
library(superheat)
library(xgboost)

setwd("~/Work/Pancreatic_Work/network_predictions")
load("assets/predicting_data.rda")

tcga <- read.table("clin_data/TCGA_mmc2.csv", sep=',', header=T, stringsAsFactors = F)
tcga$ShortBarcode <- str_sub(tcga$SampleBarcode, start=1, end=15)
rownames(tcga) <- tcga$ShortBarcode
tcga <- tcga[clin$SampleBarcode,]

c2cVar <- apply(c2cWts, 1, var, na.rm=T)

qplot(c2cVar, geom='density')
# there's a peak at 0

quantile(c2cVar[c2cVar > 0.0025], probs = 0.9)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.002501 0.007898 0.010890 0.011530 0.014460 0.048680 
#90% 
#0.01825804 


superheat(X = highVarWts, dist.method = 'manhattan', scale = F, pretty.order.rows = T, pretty.order.cols = T)




