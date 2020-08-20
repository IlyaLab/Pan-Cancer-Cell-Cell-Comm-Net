# Pan-Cancer-Cell-Cell-Comm-Net
The Pan-Cancer cell-cell communication network, data, scripts, results.

## Important files:

Data/accessories_v2.rda  -- an R data file with the following objects:
  - binMap -- bins to categorize samples by cell type score abundance
  - ccs -- the cell-cell scaffold (ccs)
  - cellECDFS -- the cell score abundance distributions

Data/CellSortedExpressonSets/cell_type_expression_tables.rda
  - The cell sorted gene expression source data that was used to create the cell-cell scaffold.
    includes: blueprintComm, blueprintRare, ramiloski, a96, a97, mab, nov, encode
  
Data/CellSortedExpressonSets/new_Celltype_map_march8.tsv
  - The map of cell names in source data to xCell data.

## Key tasks

1. Where's the data?
   - complete data is found in gs://pan-cancer-cell-cell-comm-net
   - data used to construct the scaffold: Data/CellSortedExpressionSets
     - E-MTAB-3819 (Blueprint)
     - E-MTAB-3827 (Blueprint)
     - FANTOM5 downloads hg19.gene_phase1and2combined_tpm.osc.txt.gz (https://fantom.gsc.riken.jp/5/datafiles/latest/extra/gene_level_expression/)
     - E-MTAB-3358 (FANTOM5)
     - GSE22886 (Abbas GPL96 & GPL97)
     - GSE24759 (Novershtern)
     - GSE49910 (Mabbott)
     - E-GEOD-26284 (Encode)
   - the cell-cell scaffold (network): Data/cell-cell-scaffold-narrow-march21.rda
   - xCell cell score estimates: Data/xCell_TCGA_RSEM.txt
   - Processed examples of TCGA sample data used for weighting networks (expression and xcell levels): Data/SampleExpr
   - ebpp (EB++) == GDC TCGA Pan-Cancer RNA expression data
     - RNA (Final) - EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv
     - https://gdc.cancer.gov/about-data/publications/pancanatlas

2. How was the data parsed? (it was primarily parsed to enable cluster computing)
   - probability distributions were constructed using Code/prob_model_v2.R
    -- this entails either taking a vector of cell scores over all TCGA samples, 
       and computing an empirical cumulative distribution function (ECDF)
       using the R library stats, ecdf(x) where x is the vector of cell scores.
    -- or for a given cell type, binning xcell scores into quartiles, and for 
       samples within each bin, fitting an ecdf() on expression values for
       each gene independently.
       
3. How was the cell-cell scaffold made?
   - using Code/scaffold.R
   - try it out with Example/build_scaffold_example.sh
   
4. How was the scaffold weighted for each sample?
    - using Code/compute_weights_v2.R
    - which needs, one R data file per sample with both expression and cell estimates, plus distributions & accessories (as above).
    - try a small working exmaple with Example/weight_network_example.sh
    
5. How do I get all the data?
    - The complete set of data and distributions is stored in a google bucket 
    - See Example/download_from_cloud_bucket.sh
    - Also, much of the data and results are stored in Google BigQuery tables.
    - See Analysis/SQL for examples of how to query the weighted tables and various summaries.



## Table of contents:

1. Analysis - Data and notebooks used in the analysis of networks.

2. Code - the code to build the cell-cell scaffold, weight the networks, and perform permuation testing.

3. Data - objects needed for computing the scaffold and edge weights.

4. Distributions - gene epxression and cell abundance probability distribution examples.

5. Example - examples showing how the code.

7. Results - edge weights, permuted weights, and the Differentially Weighted Edges (DWEs)



