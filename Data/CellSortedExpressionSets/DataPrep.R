library(readr)
library(stringr)
library(dplyr)


# This script takes the data sources and the map
# of cell names, and produces a table per data
# source that was then used to create the cytokine
# network scaffold. 

# The columns are individual cells and the rows
# are the bulk expression signatures.

#setwd("~/Work/Cytokine-networks")

# Cell maps
celltype_map <- read_tsv("Data/new_celltype_map_march8.tsv", col_names = FALSE)
#celltype_map$X4 <- str_replace_all(celltype_map$X3, pattern='\"', replacement = '')


# Data
xcell <- read_tsv('Data/xCell_Data/xCell_TCGA_RSEM.txt')


# Blueprint Data
blueprintRare <- read_tsv('Data/Blueprint_Data/E-MTAB-3819-query-results.tpms.tsv', comment = '#')
blueprintComm <- read_tsv('Data/Blueprint_Data/E-MTAB-3827-query-results.tpms.tsv', comment = '#')


# Fantom 5
fantom5 <- read_tsv('Data/Fantom5_Data/hg19.gene_phase1and2combined_tpm.osc.txt.gz', comment='#')
ramiloski <- read_tsv('Data/Fantom5_Data/fantom5_ExpressionLigRec.txt')
fantom5Means <- read_tsv('~/Downloads/E-MTAB-3358-query-results.tpms.tsv')
fantom5meta <-  read_tsv('~/Downloads/E-MTAB-3358-experiment-design.tsv')

# GSE22886 - GPL96
# GSE22886 - GPL97
abbas96 <- read_tsv('Data/Abbas/GSE22886-GPL96_Meta.tsv', col_names = F)
abbas97 <- read_tsv('Data/Abbas/GSE22886-GPL97_Meta.tsv', col_names = F)
abbas96Ex <- read_tsv('Data/Abbas/GSE22886-GPL96_series_matrix.txt.gz', comment = '!')
abbas97Ex <- read_tsv('Data/Abbas/GSE22886-GPL97_series_matrix.txt.gz', comment = '!')
platform96 <- read_tsv('Data/Abbas/GPL96-57554.txt', comment = '#')
platform97 <- read_tsv('Data/Abbas/GPL97-17394.txt', comment = '#')
a96 <- inner_join(platform96, abbas96Ex, by=c('ID'='ID_REF'))
a97 <- inner_join(platform97, abbas97Ex, by=c('ID'='ID_REF'))


# GSE24759
# platform
NovershternMeta <- read_tsv('Data/Novershtern/GSE24759_Meta.txt', col_names = F)
platform4685 <- read_tsv('Data/Novershtern/GPL4685-15513.txt', comment = '#')
load('Data/Novershtern/NovershternMatrix.rda')
novershtern <- cbind(data.frame(ID=rownames(NoversternMat)), as.data.frame(NoversternMat))
nov <- inner_join(platform4685, novershtern, by=c('ID'='ID'))


# GSE49910
Mabbott <- read_tsv('Data/Mabbott/GSE49910_RMA.txt.gz')
platform570 <- read_tsv('Data/Mabbott/GPL570-55999.txt', comment = '#')
mab <- inner_join(platform570, Mabbott, by=c('ID'='ID'))


# Encode
encodeMeta <- read_tsv('Data/Encode/metadata.tsv', comment = '#')
encList <- list()
# read in the array express data.
encfiles <- encodeMeta$`File accession`[encodeMeta$`Biosample term name` %in% celltype_map$X4[celltype_map$X2 == 'Encode']]
encnames <- encodeMeta$`Biosample term name`[encodeMeta$`Biosample term name` %in% celltype_map$X4[celltype_map$X2 == 'Encode']]
enctypes <- encodeMeta$`Output type`[encodeMeta$`Biosample term name` %in% celltype_map$X4[celltype_map$X2 == 'Encode']]
j <- 1
for (i in 1:length(encfiles)) {
  if (enctypes[i] == 'gene quantifications') {
    print(i)
    fileacc = encfiles[i]
    nameacc = encnames[i]
    print(nameacc)
    enctable <- read_tsv(paste0('Data/Encode/sample_level_data/', fileacc, '.tsv'))
    enctab22 <- enctable[,c('gene_id','TPM') ]
    colnames(enctab22)[2] <- nameacc
    encList[[j]] <- enctab22
    j <- j+1
  }
}

encTbl <- encList[[1]]
for (i in 2:length(encList)) {
  encTbl <- left_join(encTbl, encList[[i]])
  print(dim(encTbl))
}

encTbl$`Gene stable ID` <- unlist(lapply( X = (sapply(encTbl$gene_id, function(a) str_split(a, pattern='\\.')[1])), FUN = function(a) a[1]))

protGenes <- read_tsv('Data/protein-coding_gene.txt')
martGenes <- read_tsv('Data/mart_export.txt')

encode <- inner_join(protGenes, encTbl, by=c('ensembl_gene_id' = 'Gene stable ID'))


# unifying 'GeneSymbol'
colnames(blueprintComm)[2] <- 'GeneSymbol'
colnames(blueprintRare)[2] <- 'GeneSymbol'
colnames(ramiloski)[1] <- 'GeneSymbol'
colnames(a96)[11] <- 'GeneSymbol'
colnames(a97)[11] <- 'GeneSymbol'
colnames(mab)[11]  <- 'GeneSymbol'
colnames(nov)[9]  <- 'GeneSymbol'
colnames(encode)[2] <- 'GeneSymbol'

dataList <- list(blueprintComm, blueprintRare, fantom5, a96, a97, mab, nov, encode)

names(dataList) <- c('blueprintComm', 'blueprintRare', 'ramiloski', 'a96', 'a97', 'mab', 'nov', 'encode')

a96Colnames <- colnames(dataList[[4]])
newA96Colnames <- c()
for (i in 1:length(a96Colnames)) {
  if (a96Colnames[i] %in% abbas96[2,]) {
    idx <- which(a96Colnames[i] == abbas96[2,])
    newA96Colnames <- c(newA96Colnames, as.character(abbas96[1,idx]))
  } else {
    newA96Colnames <- c(newA96Colnames, a96Colnames[i])
  }
}
colnames(dataList[[4]]) <- newA96Colnames

a97Colnames <- colnames(dataList[[5]])
newA97Colnames <- c()
for (i in 1:length(a97Colnames)) {
  if (a97Colnames[i] %in% abbas97[2,]) {
    idx <- which(a97Colnames[i] == abbas97[2,])
    newA97Colnames <- c(newA97Colnames, as.character(abbas97[1,idx]))
  } else {
    newA97Colnames <- c(newA97Colnames, a97Colnames[i])
  }
}
colnames(dataList[[5]]) <- newA97Colnames


novColnames <- colnames(dataList[[7]])
novColnames2 <- str_split( as.character(NovershternMeta[1, ]), ',')
novColnames3 <- unlist( lapply( novColnames2, function(a) a[1]) )
newNovColnames <- c()

for (i in 1:length(novColnames)) {
  if (novColnames[i] %in% NovershternMeta[2,]) {
    idx <- which(novColnames[i] == NovershternMeta[2,])
    newNovColnames <- c(newNovColnames, as.character(NovershternMeta[1,idx]))
  } else {
    newNovColnames <- c(newNovColnames, novColnames[i])
  }
}
colnames(dataList[[7]]) <- newNovColnames

# I made the resulting tables smaller by keeping only
# cells that were in the cell map table.
save(dataList, file='cell_type_expression_tables.rda')
