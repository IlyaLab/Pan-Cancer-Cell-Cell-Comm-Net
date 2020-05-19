

# checking the scaffold

library(readr)

load('bigDataList.rda')

xcell <- read_tsv('Data/xCell_Data/xCell_TCGA_RSEM.txt')

dim(unique(scaffold))

load('Results/scaffold_median_cut.rda')

fantom5_ExpressionLigRec <- read_delim("Data/Fantom5_Data/fantom5_ExpressionLigRec.txt",
"\t", escape_double = FALSE, trim_ws = TRUE)

scaffold %>% filter(Gene == 'IL10RA')

table(scaffold$Cell, scaffold$Type)


ci <- 'Smooth muscle' # Osteoblast' # CD8+ T-cells'
print(ci)
# pull the columns
datSources <- unique(celltype_map$X2[celltype_map$X1 == ci])
datName <- as.character(unlist(sapply(datSources, function(a) dataNameMap[[a]])))
datCols <- getCols(ci, datName, dataList, celltype_map)


df = data.frame()
dataNames <- names(datCols)
for (di in dataNames) {
  print(di)
  geneVal <- datCols[[di]] %>% filter(GeneSymbol == 'IL10RA') %>% select(2)
  x <- data.frame(Source=di, Vals=datCols[[di]][,2], IL10RA=geneVal)
  colnames(x)[2] <- 'Vals'
  colnames(x)[3] <- 'IL10RA'
  df <- rbind(df, x)
}


p <- ggplot(data = df, aes(x = log(Vals), col=as.factor(Source))) + geom_density() + theme(legend.position="none")
p + geom_vline(aes(xintercept = log(IL10RA), col=as.factor(Source)))


library(RBioFabric)

#What are we missing?

xcell$X1[(!xcell$X1 %in% names(table(scaffold$Cell)))]
#"Astrocytes" "B-cells"    "CD8+ Tcm"   "CD8+ Tem"   "MPP"        "MSC"

scaffold$Cell <- as.character(scaffold$Cell)
scaffold$Cell[scaffold$Cell == 'b-cell'] <- 'B-cells'
scaffold$Cell[scaffold$Cell == 'astrocytes'] <- 'Astrocytes'

# now
#[1] "CD8+ Tcm" "CD8+ Tem" "MPP"      "MSC"

#have MPP in teh cell map
#but not the others.

