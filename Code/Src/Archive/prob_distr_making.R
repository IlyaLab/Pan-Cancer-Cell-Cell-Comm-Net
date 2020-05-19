
# building the ECDFs that are used in weighting edges

# need the network from 'netbuilding.R'
load("assets/ hybrid_cytokine_network.rda")

# list of cells that will be nodes 
load("assets/all_cells.rda")

load('step1.rda')
# need expression for each PAAD sample
load("assets/paad_expr_cell_decon_and_clin.rda")

library("AnnotationDbi")
library("org.Hs.eg.db")
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(paad), column="SYMBOL", keytype="ENTREZID", multiVals="first")

# edge genes
edges <- read.csv("network_data/LRPairs.csv",stringsAsFactors = F)
geneL <- unique(edges$Ligand.ApprovedSymbol)
geneR <- unique(edges$Receptor.ApprovedSymbol)
edgeGenes <- c(geneL, geneR)

# all samples in order
data.frame(colnames(step1)[-1], colnames(paad)[-1], clin$SampleBarcode)

# cells not present
presentCells <- allCells[!allCells %in% c("schwann", "delta")]
step1 <- step1[step1$Cell %in% presentCells,]

# create a ECDF for each cell type
library(MASS)
cellECDFS <- lapply(1:nrow(step1), function(a) {
  try( {ecdf(step1[a,-1])} )
})
names(cellECDFS) <- step1$Cell

# get some smaller matrices
# and give genes a symbol name,
step1x <- step1[,-1]
paadx <- paad[geneSymbols %in% edgeGenes,-1]
rownames(paadx) <- geneSymbols[geneSymbols %in% edgeGenes]

library(hash)

# zero function
zedf <- function(x){0}

qnames <- function(x) {
  n <- c()
  for (i in 1: (length(x)-1)) {
     n <- c(n, paste(x[i], x[i+1], sep='_'))
  }
  n
}

# top level list of cells
condECDF <- list()
cellBins <- list()

# for each cell
for (ci in presentCells) {
  
  print(ci)
  #   bin the cell quantity 
  ciquant <- quantile(x = step1x[ci,]) 

  # list with an entry for each bin
  binList <- list()
  cellBins[[ci]] <- ciquant
  
  #for each bin
  for (qi in 1:(length(ciquant)-1)) {
    print(qi)
    idx <- which(as.numeric(step1x[ci,]) > as.numeric(ciquant[qi]) & 
                   as.numeric(step1x[ci,]) <= as.numeric(ciquant[qi+1])
    )
    
    # for each gene in edgeGenes
    geneList <- list()

    for (gi in edgeGenes) {
      # fit an ECDF
      x <- na.omit(as.numeric(paadx[gi,idx]))
      if (length(x) > 3) {
        geneList[[gi]] <- ecdf(x)
      } else {
        geneList[[gi]] <- zedf
      }
    } 
    
    binList[[qi]] <- geneList
  } # end bins  
  
  names(binList) <- qnames(names(ciquant))
  condECDF[[ci]] <- binList

} # end cells



### need to know what bin to start with ###
binMap <- list()

# for each cell
for (ci in presentCells) {
  
  print(ci)
  #   bin the cell quantity 
  ciquant <- quantile(x = step1x[ci,]) 
  
  # list with an entry for cell type.
  binMap[[ci]] <- ciquant
  
}


save(paadx, step1x, presentCells, file="formated_input_paad_data.rda")

save(binMap, condECDF, cellECDFS, file="p_distr.rda")




