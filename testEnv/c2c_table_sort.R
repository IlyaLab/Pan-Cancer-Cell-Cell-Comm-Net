
# sorting the edge list

ccsdt <- data.table(ccs)

newTable <- data.table()

for (ci in allCells) {
  for (cj in allCells) {
    x <- subset(ccsdt, cellL == ci & cellR == cj)
    newTable <- rbind(newTable, x)
  }
}

