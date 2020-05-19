

# weighted graphs for each sample #
library(readr)
library(data.table)
library(foreach)

#load('Results/cell-cell-scaffold-narrow-march11.rda')  ## ccs: cell to cell scaffold

load('Data/ebpp.rda')  ## PanCancer expression data

xcell <- read_tsv('Data/xCell_Data/xCell_TCGA_RSEM.txt')  ## the cell scores

load('Results/p_distr.rda')  ## the probability model

# preproc

#ccs <- apply(ccs, 2, as.character)
xcell <- data.table(xcell)
ebpp <- cbind(data.frame(GeneID=rownames(ebpp), stringsAsFactors=F), ebpp)
ebpp <- data.table(ebpp)

colnames(xcell) <- str_replace_all(colnames(xcell), pattern='\\.', replacement='-')

xcell$X1[3] <- 'astrocytes'
xcell$X1[4] <- 'b-cell'

gc() 

# > names(condECDF[[ci]])
# [1] "0%_25%"   "25%_50%"  "50%_75%"  "75%_100%"

binfun <- function(x, ci, binMap) {
  # given a cell quantity for a given sample,
  # what bin was it placed into?
  # a <- min(which(binMap[[ci]] < x))
  # (a-1)
  if (x >= last(binMap[[ci]])) {
    a <- length(binMap[[ci]])-1
  } else if (x <= first(binMap[[ci]])) {
    a <- 1
  } else {
    a <- max(which(binMap[[ci]] < x))
  }
  return(a)
}


res0 <- data.frame(row.names = 1:nrow(ccs))


  # for each edge in the cell-to-cell network
for (si in 2:ncol(xcell)) {

  barcode <- str_sub(colnames(ebpp)[si], start=1, end=15)

  pList <- list()

  xcell1 <- xcell[,barcode,with=F]
  ebpp1  <- ebpp[,si,with=F]

  for (c1 in unique(xcell$X1)) {
      # get the look up values for each element in the prob. statement
      c1val <- xcell[X1==c1, barcode, with=F]
      bin1 <- binfun(c1val, c1, binMap)  # get the bin for this cell type
      p_l1 <- condECDF[[c1]][[bin1]]     # get the list of ECDFs
      pList[[c1]] <- p_l1                # write it to the list
      # a sample will always be in the same bin for a given cell type
  } # end cell loop

  # write to
  # /titan/cancerregulome9/workspaces/users/dgibbs/CytokineNetworkData
  pathName <- paste0('/titan/cancerregulome9/workspaces/users/dgibbs/CytokineNetworkData/',
                    'sample_',si,'.rda')
  save(xcell1, ebpp1, pList, file=pathName)

}



