library(dplyr)

# use the LM22 network,  add in single cell pancreas edges
# because that's how we deconvoluted.

# import data
load("~/Work/Cytokine_Scaffold/cibersort_network/data/lm22.full.063017.RData")
load("single_cell_data/pancreatic_cells.rda")
colnames(h) <- h_ann$cell_type1

# cells
cellsNames <- unique(colnames(lm22.full))
pancCells <- c("acinar","activated_stellate","alpha","beta","delta","ductal","endothelial","epsilon","gamma","quiescent_stellate","schwann")

# summarize panc cells
pancDF <- data.frame(row.names = rownames(h))
for (pci in pancCells) {
  print(pci)
  hci <- h[,colnames(h) == pci]
  cellMean <- apply(hci, 1, function(a) mean(a[a > 0], na.rm=T))
  pancDF <- cbind(pancDF, data.frame(pci = cellMean))
}
colnames(pancDF) <- pancCells

# #############################################################
# Now, we need to add in the extra edges
# and then we build the regular LM22 network
# but we add a new loop of building the pancreas cells in.
# will need a different threshold.
# then we will need the new edge weight algorithm.


# edge genes
edges <- read.csv("network_data/LRPairs.csv",stringsAsFactors = F)
geneL <- unique(edges$Ligand.ApprovedSymbol)
geneR <- unique(edges$Receptor.ApprovedSymbol)
edgeGenes <- c(geneL, geneR)

# subset the gene expression data
sublm22 <- lm22.full[rownames(lm22.full) %in% edgeGenes,]

# missing some
geneL <- geneL[geneL %in% rownames(sublm22)]
geneR <- geneR[geneR %in% rownames(sublm22)]
edges <- edges[edges$Ligand.ApprovedSymbol %in% geneL & edges$Receptor.ApprovedSymbol %in% geneR, ]

immThresh <- 250
pancThresh <- 5

resnet <- data.frame(stringsAsFactors = F)
for (ci in cellsNames) {
  #
  print(ci)
  for (li in geneL) {
    cexpr <- sublm22[li,colnames(sublm22) == ci] 
    if (any(!is.na(cexpr)) & mean(cexpr, na.rm = T) > immThresh) {
      resnet <- rbind(resnet, data.frame(From=ci, To=li, Type='L', stringsAsFactors = F))
    }
  }
  for (ri in geneR) {
    cexpr <- sublm22[ri,colnames(sublm22) == ci] 
    if (any(!is.na(cexpr)) & mean(cexpr, na.rm = T) > immThresh) {
      resnet <- rbind(resnet, data.frame(From=ri, To=ci, Type='R', stringsAsFactors = F))
    }
  }
}

# add panc cells here # 
for (ci in pancCells) {
  #
  print(ci)
  for (li in geneL) {
    cexpr <- pancDF[li,colnames(pancDF) == ci]
    if (any(!is.na(cexpr)) & mean(cexpr, na.rm = T) > pancThresh) {
      resnet <- rbind(resnet, data.frame(From=ci, To=li, Type='L', stringsAsFactors = F))
    }
  }
  for (ri in geneR) {
    cexpr <- pancDF[ri,colnames(pancDF) == ci] 
    if (any(!is.na(cexpr)) & mean(cexpr, na.rm = T) > pancThresh) {
      resnet <- rbind(resnet, data.frame(From=ri, To=ci, Type='R', stringsAsFactors = F))
    }
  }
}


for (i in 1:nrow(edges)) {
  li <- edges[i,"Ligand.ApprovedSymbol"]
  ri <- edges[i,"Receptor.ApprovedSymbol"]
  resnet <- rbind(resnet, data.frame(From=li, To=ri, Type='M', stringsAsFactors = F))
}

hybridNet <- resnet
save(hybridNet, file="hybrid_cytokine_network.rda")
write.table(hybridNet, file="hybrid_cytokine_network.tsv", sep='\t', row.names = F, quote=F)
