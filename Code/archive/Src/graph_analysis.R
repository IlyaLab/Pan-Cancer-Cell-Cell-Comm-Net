

# graph analysis


library(igraph)
library(RBioFabric)

load('scaffold_median_cut.rda')
edges <- read_tsv('Data/Fantom5_Data/fantom5_PairsLigRec.txt')

g1 <- igraph::graph_from_edgelist(matrix(data=c(scaffold$Cell, scaffold$Gene), byrow = F, ncol=2), directed = T)

g2 <- igraph::graph_from_edgelist(matrix(data=c(edges$Ligand.ApprovedSymbol, edges$Receptor.ApprovedSymbol), ncol=2, byrow=F), directed=T)

g <- g1+g2
plot(g)

RBioFabric::bioFabric(g)

plot.igraph(g, vertex.size=0.5, vertex.label=F, edge.arrow.size=0.5)
