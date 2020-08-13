
# building the ECDFs that are used in weighting edges
library(hash)
library(MASS)
library(readr)
library(parallel)
library(stringr)

# weighted graphs for each sample #

load('Data/cell-cell-scaffold-narrow-march21.rda')

load('Data/ebpp.rda')

xcell <- read_tsv('Data/xCell_Data/xCell_TCGA_RSEM.txt')

ccs$cellL <- as.character(ccs$cellL)
ccs$cellR <- as.character(ccs$cellR)
ccs$ligand <- as.character(ccs$ligand)
ccs$receptor <- as.character(ccs$receptor)

xcell$X1[xcell$X1 == "B-cells"] <- 'b-cell'
xcell$X1[xcell$X1 == "Astrocytes"] <- 'astrocytes'


# list of cells that will be nodes
allCells <- unique(xcell$X1)
all(allCells %in% c(ccs$cellL, ccs$cellR))
#true

edgeGenes <- unique(c(ccs$ligand, ccs$receptor))

ebpp2 <- ebpp[edgeGenes,]

rm(ebpp)
gc()

# create a ECDF for each cell type
cellECDFS <- lapply(1:nrow(xcell), function(a) {
  print(xcell$X1[a]);
  try( {ecdf(as.numeric(xcell[a,-1]))} )
})
names(cellECDFS) <- xcell$X1


# zero function
zedf <- function(x){0}

qnames <- function(x) {
  n <- c()
  for (i in 1: (length(x)-1)) {
     n <- c(n, paste(x[i], x[i+1], sep='_'))
  }
  n
}


### need to know what bin to start with ###
binMap <- list()

# for each cell
for (ci in allCells) {

  print(ci)
  #   bin the cell quantity
  vals <- as.numeric(xcell[xcell$X1 == ci,-1])
  ciquant <- quantile(x = vals)

  # list with an entry for cell type.
  binMap[[ci]] <- ciquant
}


# top level list of cells
condECDF <- list()
cellBins <- list()

# EB++ barcodes starts at 1
shortBarcodes <- str_sub(colnames(ebpp2), 1, 15)
# xcell barcodes starts at 2
xcellBarcodes <- colnames(xcell)[-1]
xcellBarcodes <- str_replace_all(xcellBarcodes, '\\.', '-')

# for each cell
for (ci in allCells) {

  print(ci)
  #   bin the cell quantity
  vals <- as.numeric(xcell[xcell$X1 == ci,-1])
  ciquant <- quantile(x = vals)

  # list with an entry for each bin
  binList <- list()
  cellBins[[ci]] <- ciquant

  #for each bin
  for (qi in 1:(length(ciquant)-1)) {
    print(qi)

    # get the samples in this quantile
    idx <- which(vals > as.numeric(ciquant[qi]) &   # ## !! THIS IS INDEX TO xCELL !
                 vals <= as.numeric(ciquant[qi+1]))

    # and get the barcodes and index to expression data
    qiXcellBarcodes <- xcellBarcodes[idx]
    jdx <- which(shortBarcodes %in% qiXcellBarcodes)  # get index for ebpp

    # for each gene in edgeGenes
    geneList <- mclapply(X=edgeGenes, FUN = function(gi) {
      # get the values for this quantile
      x <- na.omit(as.numeric(ebpp2[gi,jdx]))
      # and make a ecdf
      if (length(x) > 3) {
        ecdf(x)
      } else {
        zedf
      }
    }, mc.silent = T, mc.cores = 5, mc.preschedule = T)

    names(geneList) <- edgeGenes

    binList[[qi]] <- geneList
  } # end bins

  names(binList) <- qnames(names(ciquant))
  condECDF[[ci]] <- binList

} # end cells



save(binMap, condECDF, cellECDFS, file="p_distr_v2.rda")

# for each cell
#    for each quantile
#       for each gene
