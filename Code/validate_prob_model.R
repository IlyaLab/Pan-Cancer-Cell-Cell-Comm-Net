

# weighted graphs for each sample #

load('p_distr.rda')

### validation ###
allCells <- names(cellECDFS)
allGenes <- names(condECDF[[1]][[1]])

for (ci in allCells) {

  ad <- cellECDFS[[ci]]
  bd <- condECDF[[ci]]

  if (length(bd) < 4) {
    print(paste0("error, condECDF didn't have 4 quantiles:  ", ci))
  }

  if ( max(summary(ad)) <= 0) {
    print(paste0('error, cellECDF has max zero:  ', ci))
  }

  for (j in 1:4) {
    cd <- bd[[j]]

    zeroCnt <- 0
    for (g in 1:length(cd)) {

      dd <- cd[[g]]

      if (any(class(dd) == 'ecdf')) {  # could have a zero ecdf fun

        if ( max(summary(dd)) <= 0 ) {
          print(paste0('error, cellECDF has max zero:  ', ci, ' ', j, ' ', g))
        }

      } else {
        zeroCnt <- zeroCnt + 1  # check if we have too many of these
      }

    }
    if (zeroCnt > 100) {
      print(paste0('error, zero count > 100:  ', ci, ' ', j))
    }
  }
}

print('done')

