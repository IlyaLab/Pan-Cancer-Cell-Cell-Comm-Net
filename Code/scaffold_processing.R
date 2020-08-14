
library(readr)
library(dplyr)

# building full network
load('scaffold_median_cut_aug13.rda')

# edge genes
edges <- read_tsv("Data/fantom5_PairsLigRec.txt")


### Making the full cell-L-R-cell graph ####

edgeTbl <- edges %>% select(Ligand.ApprovedSymbol, Receptor.ApprovedSymbol) %>% rename(Ligand=Ligand.ApprovedSymbol, Receptor=Receptor.ApprovedSymbol)

a <- scaffold %>% filter(Type == 'Ligand')
b <- scaffold %>% filter(Type == 'Receptor')

leftSide <- left_join(x=a, y=edgeTbl, by=c('Gene'='Ligand'))
ccScaffold <- left_join(x=leftSide, y=b, by=c('Receptor' = 'Gene'))

# a NA cell was found in the CellR column
ccScaffold <- na.omit(ccScaffold)

save(ccScaffold, file='cell-cell-scaffold_aug13.rda')

### Just the most used columns ###

ccs <- ccScaffold %>% select(Cell.x, Gene, Receptor, Cell.y) %>% rename(cellL=Cell.x, ligand=Gene, receptor=Receptor, cellR=Cell.y)
save(ccs, file='cell-cell-scaffold-narrow-aug13.rda')


### Making a long graph representation ###

ccs$cellL <- as.character(ccs$cellL)
ccs$cellR <- as.character(ccs$cellR)
ccs$ligand <- as.character(ccs$ligand)
ccs$receptor <- as.character(ccs$receptor)

ccsLong <- unique(data.frame(From = c(ccs$cellL, ccs$ligand, ccs$receptor), To=c(ccs$ligand, ccs$receptor, ccs$cellR), stringsAsFactors = F))
ccsSif <- cbind(From=ccsLong$From, Wt=1, To=ccsLong$To)
write.table(ccsSif, file='ccsSif_aug13.sif', sep='\t', row.names = F, col.names = F, quote=F)


### Making a two part graph C-L  ->>  R-C  ###

part1 <- sapply(1:nrow(ccs), function(i) paste0(ccs$cellL[i], '_', ccs$ligand[i]))
part2 <- sapply(1:nrow(ccs), function(i) paste0(ccs$receptor[i], '_', ccs$cellR[i]))
cct <- data.frame(From=part1, To=part2, stringsAsFactors = F)

