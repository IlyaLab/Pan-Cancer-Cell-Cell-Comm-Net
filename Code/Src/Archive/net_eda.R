
load("all_cells.rda")
load("hybrid_cytokine_network.rda")

allCells <- c(pancCells, cellsNames)

outNet <- hybridNet[hybridNet$From %in% allCells,]

table(outNet$From)


inNet <- hybridNet[hybridNet$To %in% allCells,]

table(inNet$To)

