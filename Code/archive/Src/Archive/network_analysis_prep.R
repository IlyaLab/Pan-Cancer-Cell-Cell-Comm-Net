#cell 2 cell analysis

library(dplyr)
library(tidyr)
library(ggridges)
library(stringr)

clin <- read.table("clin_data/pancan_clinical-20180322-124211_no_mets.csv", sep=',', stringsAsFactors = F, header=T)
resp <- read.table("clin_data/panimmune_clinical-20180322-115051.csv", sep=',', stringsAsFactors = F, header = T)
load("assets/c2c_wts.rda")

ids <- clin$bcr_patient_barcode
ids <- str_join(ids, '-01', sep='')
clin$SampleBarcode <- ids

subids <- intersect(clin$SampleBarcode, colnames(res0))

ids <- resp$ParticipantBarcode
ids <- str_join(ids, '-01', sep='')
resp$SampleBarcode <- ids

res1 <- res0[,colnames(res0) %in% subids]
clin <- clin[clin$SampleBarcode %in% subids,]
resp <- resp[resp$SampleBarcode %in% subids,]

idx <- match(table = clin$SampleBarcode, x = colnames(res1))
c2 <- clin[idx,]
all(c2$SampleBarcode == colnames(res1))

idx <- match(table = resp$SampleBarcode, x = colnames(res1))
r2 <- resp[idx,]
all(r2$SampleBarcode == colnames(res1))

clin <- c2
resp <- r2
rm(c2, r2, ids, subids, idx)

paths <- lapply(1:nrow(c2cNet), function(i) paste(c2cNet[i,1], c2cNet[i,2], c2cNet[i,3], c2cNet[i,4], sep='_'))
c2cNet <- cbind(c2cNet, data.frame(Path=unlist(paths)))

res2 <- cbind(data.frame(Path=unlist(paths)), data.frame(res1))
resLong <- gather(data = res2, key=Barcodes, value = Score, -Path)

c2cWts <- res1
c2cLong <- resLong

rm(res2, paths, res1, resLong)

save(list=ls(), file="assets/predicting_data.rda")

