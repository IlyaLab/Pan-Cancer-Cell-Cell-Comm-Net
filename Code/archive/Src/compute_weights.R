# weighted graphs for each sample #
library(hash)
library(data.table)
library(foreach)

workFun <- function(si) {
  # for each edge in the cell-to-cell network
  # hold a value for each edge
  # load the expression and xcell for this sample
  filepath = '/titan/cancerregulome9/workspaces/users/dgibbs/CytokineNetworkData/'
  load(paste0(filepath, 'cell-cell-scaffold-narrow-march21.rda'))  ## ccs: cell to cell scaffold
  load(paste0(filepath, 'accessories.rda'))  ## ccs: cell to cell scaffold

  vals <- numeric(nrow(ccs))

  try({

    curr1 <- ''
    curr2 <- ''

    load(paste0(filepath, 'Expr/sample_',si ,'.rda')) # expr


    for (ei in 1:nrow(ccs)) {  

	try ({

        # get the look up values for each element in the prob. statement
        c1 <- as.character(ccs$cellL[ei]) # cell 1
        c2 <- as.character(ccs$cellR[ei]) # cell 2
        li <- as.character(ccs$ligand[ei]) # ligand
        re <- as.character(ccs$receptor[ei]) # receptor

        c1val <- xcell1[cellIdx[[c1]], 2]
        c2val <- xcell1[cellIdx[[c2]], 2]

        lival <- ebpp1[geneIdx[[li]], 2]
        reval <- ebpp1[geneIdx[[re]], 2]

        bin1 <- cellBin1[[c1]] # look up bin here
        bin2 <- cellBin1[[c2]] # look up bin here

        print(paste(c1,c2, li, re, c1val, c2val, lival, reval, bin1, bin2, sep=' '))

        if (curr1 != paste0(c1,bin1)) {
          # then load it.
          load(paste0(filepath, 'Prob/prob_',c1,'_',bin1,'.rda'))
          p1 <- x
          curr1 <- paste0(c1,bin1)
        }
        
        if (curr2 != paste0(c2,bin2)) {
          # then load it.
          load(paste0(filepath, 'Prob/prob_',c2,'_',bin2,'.rda'))
          p2 <- x
          curr2 <- paste0(c2,bin2)
        }

        # use a hash to look up these
        p_l1 <- p1[[li]](lival)
        p_c1 <- cellECDFS[[c1]](c1val)
        p_r2 <- p2[[re]](reval)
        p_c2 <- cellECDFS[[c2]](c2val)

        # save the result of the product
        vals[ei] <- p_l1*p_c1*p_r2*p_c2

	})# end edge try

    } # end edge loop

  }) # end try

  return(vals)
}

## 11069 total samples


cl <- parallel::makeForkCluster(20)

for (si in seq(from=1, to=400, by=200)) { # to=11069, by=200)) {

  sampls  <- seq(from=si, to=(si+199))

  doParallel::registerDoParallel(cl)
  parRes <- foreach(i = sampls) %dopar% {
    workFun(i)
  }

  save(sampls, parRes, file=paste0("Results/parRes_",si,".rda"))
  parRes <- list()
  gc()

}

