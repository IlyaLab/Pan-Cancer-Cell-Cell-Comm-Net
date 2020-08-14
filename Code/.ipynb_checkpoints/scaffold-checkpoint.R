
library(readr)
library(stringr)
library(ggplot2)
library(dplyr)

# Fantom5 edges
edges <- read_tsv('Data/fantom5_PairsLigRec.txt')

# Cell maps
celltype_map <- read_tsv("Data/CellSortedExpressionSets/new_Celltype_map_march8.tsv")

# Cell data
xcell <- read_tsv('Data/xCell_TCGA_RSEM.txt')

#load('bigDataList.rda')
load('Data/CellSortedExpressionSets/cell_type_expression_tables.rda') # dataList2

#                      1              2              3      4    5    6    7     8
#dataList <- list(blueprintComm, blueprintRare, ramiloski, a96, a97, mab, nov, encode)


table(celltype_map$X2)
celltype_map$X2[celltype_map$X2 == 'Fantom 5'] <- 'Fantom5'
celltype_map$X2[celltype_map$X2 == 'Blueprint'] <- 'blueprint'
celltype_map$X2[celltype_map$X2 == 'GSE22886 Abbas'] <- 'Abbas'

dataNameMap <- list()
dataNameMap[['blueprint']] <- c('blueprintComm', 'blueprintRare')
dataNameMap[['Fantom5']] <- c('ramiloski')
dataNameMap[['Abbas']] <- c('a96', 'a97')
dataNameMap[['Mabbott']] <- c('mab')
dataNameMap[['Noverstern']] <- c('nov')
dataNameMap[['Encode']] <- c('encode')

condensedCellMap <- unique(celltype_map[,c(1,2,4)])

# fix nov col names
novCols <- colnames(dataList2[[7]])
novCols2 <- str_split(novCols, ',')
novCols3 <- lapply(novCols2, function(a) a[1])
novCols <- unlist(novCols3)
colnames(dataList2[[7]]) <- novCols
rm(novCols2, novCols3)

colnames(dataList2[['encode']])[2] <- 'GeneSymbol'
colnames(dataList2[[3]])[1] <- 'GeneSymbol'

get10p <- function(x, pi) {
  # return 10% of the range.
  #maxx <- max(x[,2], na.rm = T)
  #minx <- min(x[,2], na.rm = T)
  #return(minx + (0.1 * (maxx-minx)))
  colnames(x)[2] <- 'vals'
  y <- x %>% dplyr::filter(vals > 0) %>% select(vals) %>% first()
  quantile(y, probs = pi, na.rm = T)
}


getCols <- function(celli, src, dList, cmap) {
  res0 <- list()
  cellnames <- cmap[cmap$X1 == ci, 4]$X4 # get the right cell type names
  for (si in src) { # for each source
    print(paste0('     ',si))
    idx <- which(colnames(dList[[si]]) %in% cellnames) # get the columns for this cell ID
    if (length(idx) > 0) {
      for (i in idx) {
        thiscellname <- colnames(dList[[si]])[i]
        x <- dList[[si]][,c('GeneSymbol', thiscellname)]  # get each col with the gene symbol
        res0[[paste0(si,'_',celli,'_',i)]] <- x
      }
    }
  }
  res0
}


makeCut <- function (i, datMats, cuts, geneSym) {
    if (geneSym %in% datMats[[i]]$GeneSymbol) {
      a <- datMats[[i]]
      b <- a %>% filter(GeneSymbol == geneSym) %>% select(2) %>% first()
      if (any(sapply(b, function(aaa) {!is.na(aaa)}))) {
        bavg <- mean(b, na.rm = T)
        if (bavg > cuts[[i]][1]) {
          return(sum(bavg > a[,2], na.rm = T) / nrow(a)) # true in upper range
        } else {
          return(FALSE)
        }
      }
    } else {
      return(NA)
  }
    return(-999) # ERROR
}


# WANT:  Cell - Molecule pairs.
# is it expressed?

scaffold <- data.frame()

# for each cell type:
for (ci in unique(celltype_map$X1)) {
  print(ci)
  # pull the columns
  datSources <- unique(celltype_map$X2[celltype_map$X1 == ci])
  datName <- as.character(unlist(sapply(datSources, function(a) dataNameMap[[a]])))
  datCols <- getCols(ci, datName, dataList2, celltype_map)

  # get the quantile to cut
  cutLevel <- 0.5
  datCuts <- lapply(datCols, function(a) get10p(a, cutLevel))

  # for each gene
  for (ei in 1:nrow(edges)) {
    try({
      genes <- c(edges$Ligand.ApprovedSymbol[ei], edges$Receptor.ApprovedSymbol[ei])

      for (gi in 1:2) {
        # we get the number of samples where this gene is expressed above the 10% cutoff.
        vals <- unlist(lapply(1:length(datCols), function(i) makeCut(i, datCols, datCuts, genes[gi])))

        if ( (sum(vals > 0, na.rm = T)/length(vals)) > 0.5) {  # if more than one sample shows expression of this molecule.
          if (gi == 1) { # then it's a ligand
            scaffold <- rbind(scaffold, data.frame(Cell=ci, Gene=genes[gi], Type='Ligand', Wt=median(vals, na.rm = T), Vote=sum(vals > 0)/length(vals), Samples=length(vals)))
          } else {
            scaffold <- rbind(scaffold, data.frame(Cell=ci, Gene=genes[gi], Type='Receptor', Wt=median(vals, na.rm = T), Vote=sum(vals > 0)/length(vals), Samples=length(vals)))
          }
        }
      } # end genes
    })
  } # end edges
} # end cells


scaffold <- unique(scaffold)
save(scaffold, file='scaffold_median_cut_aug13.rda')


