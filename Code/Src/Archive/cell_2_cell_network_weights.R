# weighted networks for each sample

# condECDF
# cell -> bin -> gene -> ECDF

# cellECDF
# cell -> ECDF

# want 
# P(L | C) p(C)  * P(R | C) P(C)
# p_l1    p_c1     p_r2     p_c2


library(dplyr)
library(parallel)

# need the network
load("hybrid_cytokine_network.rda")
load("p_distr.rda")
load("formated_input_paad_data.rda")

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

a_net <- hybridNet[hybridNet$Type == 'L', c(1,2)]
b_net <- hybridNet[hybridNet$Type == 'M', c(1,2)]
c_net <- hybridNet[hybridNet$Type == 'R', c(1,2)]

d_net <- inner_join(a_net, b_net, by =c('To' = 'From'))
e_net <- inner_join(b_net, c_net, by =c('To' = 'From'))
colnames(d_net) <- c('From','Lig','Rec')
colnames(e_net) <- c('Lig','Rec','To')

# cell to cell net
c2cNet <- inner_join(d_net, e_net)

rm(a_net, b_net, c_net, d_net, e_net)


res0 <- data.frame(row.names = 1:nrow(c2cNet))

workFun <- function(si, c2cNet, step1x, paadx, binfun, binMap, condECDF, cellECDF) {
  # for each edge in the cell-to-cell network
  
  vals <- replicate(expr = 0, n = nrow(c2cNet))
  
  for (ei in 1:nrow(c2cNet)) {
    if (ei %% 200 == 0) {print(ei)}
    # get the look up values for each element in the prob. statement
    c1 <- c2cNet[ei, 1]
    c2 <- c2cNet[ei, 4]
    li <- c2cNet[ei, 2]
    re <- c2cNet[ei, 3]
    c1val <- step1x[c1, si]
    c2val <- step1x[c2, si]
    lival <- paadx[li, si]
    reval <- paadx[re, si]
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
  }
  return(vals)  
}

# for each sample
resP <- parallel::mclapply(1:ncol(paadx), FUN = function(i) {workFun(i, c2cNet, step1x, paadx, binfun, binMap, condECDF, cellECDF)}, mc.cores = 16)

save(resP, file="resP_results.rda")

res0 <- do.call("cbind", resP)


colnames(res0) <- colnames(paadx)
save(res0, file="c2c_wts.rda")


