

# weighted graphs for each sample #
library(data.table)
library(parallel)


# > names(condECDF[[ci]])
# [1] "0%_25%"   "25%_50%"  "50%_75%"  "75%_100%"


workFun <- function(si, ccs, xcell, ebpp, binfun, binMap, condECDF, cellECDF) {
  # for each edge in the cell-to-cell network
  curr1 <- ''
  curr2 <- ''

  # load the expression and xcell for this sample
  filepath = '/titan/cancerregulome9/workspaces/users/dgibbs/CytokineNetworkData/'
  load(paste0(filepath, 'Expr/sample_',si ,'.rda')) # expr
  load(paste0(filepath, 'cell-cell-scaffold-narrow-march11.rda'))  ## ccs: cell to cell scaffold
  load(paste0(filepath, 'accessories.rda'))  ## ccs: cell to cell scaffold

  # hold a value for each edge
  vals <- numeric(nrow(ccs))

  for (ei in 1:100) { # nrow(ccs)) {

    try({
      # get the look up values for each element in the prob. statement
      c1 <- as.character(ccs[ei, 1]) # cell 1
      c2 <- as.character(ccs[ei, 4]) # cell 2
      li <- as.character(ccs[ei, 2]) # ligand
      re <- as.character(ccs[ei, 3]) # receptor

      c1val <- xcell1[xcell1$X1==c1, 2]
      c2val <- xcell1[xcell1$X1==c2, 2]

      lival <- ebpp1[ebpp1$GeneID==li, 2]
      reval <- ebpp1[ebpp1$GeneID==re, 2]

      bin1 <- cellBin1[[c1]] # look up bin here
      bin2 <- cellBin1[[c2]] # look up bin here

      if (curr1 != paste0(c1,bin1)) {
        # then load it.
        load(paste0(filepath, 'Prob/prob_',c1,'_',bin1,'.rda'))
        p1 <- x
      }
      if (curr2 != paste0(c2,bin2)) {
        # then load it.
        load(paste0(filepath, 'Prob/prob_',c1,'_',bin1,'.rda'))
        p2 <- x
      }

      p_l1 <- p1[[li]](lival)
      p_c1 <- cellECDFS[[c1]](c1val)
      p_r2 <- p2[[re]](reval)
      p_c2 <- cellECDFS[[c2]](c2val)
      # save the result of the product
      vals[ei] <- p_l1*p_c1*p_r2*p_c2
    }) # end try

  } # end edge loop

  return(vals)

}


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


