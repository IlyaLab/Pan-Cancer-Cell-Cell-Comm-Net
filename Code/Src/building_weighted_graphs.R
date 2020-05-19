

# weighted graphs for each sample #
library(readr)
library(data.table)
library(parallel)

load('Results/cell-cell-scaffold-narrow-march11.rda')  ## ccs: cell to cell scaffold

load('Data/ebpp.rda')  ## PanCancer expression data

xcell <- read_tsv('Data/xCell_Data/xCell_TCGA_RSEM.txt')  ## the cell scores

load('Results/p_distr.rda')  ## the probability model

# preproc

ccs <- apply(ccs, 2, as.character)
xcell <- data.table(xcell)
ebpp <- cbind(data.frame(GeneID=rownames(ebpp), stringsAsFactors=F), ebpp)
ebpp <- data.table(ebpp)

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

workFun <- function(si, ccs, xcell, ebpp, binfun, binMap, condECDF, cellECDF) {
  # for each edge in the cell-to-cell network


  vals <- lapply(1:nrow(ccs), function(ei) {
    try({
      # get the look up values for each element in the prob. statement
      c1 <- ccs[ei, 1] # cell 1
      c2 <- ccs[ei, 4] # cell 2
      li <- ccs[ei, 2] # ligand
      re <- ccs[ei, 3] # receptor

      c1val <- xcell[X1==c1, si, with=F]
      c2val <- xcell[X1==c2, si, with=F]

      lival <- ebpp[GeneID==li, si, with=F]
      reval <- ebpp[GeneID==re, si, with=F]

      if (all(!sapply(c(c1val, c2val, lival, reval), is.na))) {
        bin1 <- binfun(c1val, c1, binMap)
        bin2 <- binfun(c2val, c2, binMap)
        p_l1 <- condECDF[[c1]][[bin1]][[li]](lival)
        p_c1 <- cellECDFS[[c1]](c1val)
        p_r2 <- condECDF[[c2]][[bin2]][[re]](reval)
        p_c2 <- cellECDFS[[c2]](c2val)
        # save the result of the product
        vals[ei] <- p_l1*p_c1*p_r2*p_c2
      }
    }) # end try
  }) # end lapply
  return(vals)

}


ccsSub <- ccs[1:1000,]

cl <- parallel::makeForkCluster(4)
doParallel::registerDoParallel(cl)
parRes <- foreach(i = 2:20, .combine = 'c') %dopar% {
  workFun(i, ccsSub, xcell, ebpp, binfun, binMap, condECDF, cellECDF)
}

# for each sample
#resP <- parallel::mclapply(1:ncol(ebpp), FUN = function(i) {workFun(i, c2cNet, step1x, paadx, binfun, binMap, condECDF, cellECDF)}, mc.cores = 10)
#resP <- parallel::mclapply(2:ncol(ebppSubset), FUN = function(i) {workFun(i, ccs, xcell, ebppSubset, binfun, binMap, condECDF, cellECDF)}, mc.cores = 3)

save(parRes, file="weighted_c2c_graph_list.rda")

#res0 <- do.call("cbind", resP)

#colnames(res0) <- colnames(ebpp)
#save(res0, file="c2c_wts.rda")


