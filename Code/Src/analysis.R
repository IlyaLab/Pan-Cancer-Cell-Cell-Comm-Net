### analysis ###

> head(ccs)
cellL ligand receptor        cellR
1   aDC    A2M     LRP1   astrocytes
2   aDC    A2M     LRP1       b-cell
3   aDC    A2M     LRP1     CD4+ Tem
4   aDC    A2M     LRP1 CD8+ T-cells
5   aDC    A2M     LRP1          cDC
6   aDC    A2M     LRP1 Chondrocytes
> sum(unique(ccs$ligand) %in% rownames(ebpp))
[1] 603
> length(unique(ccs$ligand))
[1] 607
> sum(unique(ccs$receptor) %in% rownames(ebpp))
[1] 650
> length(unique(ccs$receptor))
[1] 661

> unique(ccs$receptor)[!(unique(ccs$receptor) %in% rownames(ebpp))]
[1] "ACKR2"  "ACKR4"  "ACKR3"  "SLC9C2" "NPY4R"  "C5AR2"  "PTGDR2" "GP1BB"  "IFNLR1" "IGFLR1" "ASIC3"

unique(ccs$ligand)[!(unique(ccs$ligand) %in% rownames(ebpp))]
[1] "PKM"     "IFNL1"   "C1QTNF5" "UTS2B"

