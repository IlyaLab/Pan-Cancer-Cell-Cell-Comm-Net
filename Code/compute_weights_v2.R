#!/usr/bin/env Rscript

# weighted graphs for each sample #
library(hash)
library(data.table)

workFun <- function(si, filepath) {
  # for each edge in the cell-to-cell network
  # hold a value for each edge
  # load the expression and xcell for this sample
  filepath = '/titan/cancerregulome9/workspaces/users/dgibbs/CytokineNetworkData/'
  # filepath = '/home/davidgibbs/Work/Cytokine-networks/compute/'

  load(paste0(filepath, 'accessories_v2.rda'))  ## ccs: cell to cell scaffold

  vals <- numeric(nrow(ccs))

  curr1 <- ''
  curr2 <- ''

  try ({

  #load(paste0(filepath, 'Expr_v2/sample_',si ,'.rda')) # expr
  load(si)

  barcode <- colnames(ebpp1)[2]
  ebpp1 <- as.matrix(ebpp1[,2])
  xcell1 <- as.matrix(xcell1[,2])

  for (ei in 1:nrow(ccs)) {

	try ({

        # get the look up values for each element in the prob. statement
        c1 <- ccs[ei,1] # cell 1
        li <- ccs[ei,2] # ligand
        re <- ccs[ei,3] # receptor
        c2 <- ccs[ei,4] # cell 2

        c1val <- xcell1[cellIdx[[c1]], 1]
        c2val <- xcell1[cellIdx[[c2]], 1]

        lival <- ebpp1[geneIdx1[[li]], 1] # sample specific geneIdx1
        reval <- ebpp1[geneIdx1[[re]], 1]

        bin1 <- cellBin1[[c1]] # look up bin here
        bin2 <- cellBin1[[c2]] # look up bin here

        if (curr1 != paste0(c1,bin1)) {
          # then load it.
          load(paste0(filepath, 'Prob_v2/prob_',c1,'_',bin1,'.rda'))
          p1 <- x
          curr1 <- paste0(c1,bin1)
        }

        if (curr2 != paste0(c2,bin2)) {
          # then load it.
          load(paste0(filepath, 'Prob_v2/prob_',c2,'_',bin2,'.rda'))
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

	}) # end try

  } # end edge loop

  }) # end try


  res0 <- list(Barcode=barcode, Vals=vals)
  return(res0)
}


args = commandArgs(trailingOnly=TRUE)

filepath = '/titan/cancerregulome9/workspaces/users/dgibbs/CytokineNetworkData/'

si <- args[1]  ## file name of data

ti <- args[2]

res1 <- workFun(si, filepath)

fileout <- paste0(filepath,"Results_v2/res_",ti,".rda")

save(si, ti, res1, file=fileout)

