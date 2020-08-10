# Pan-Cancer-Cell-Cell-Comm-Net
The Pan-Cancer cell-cell communication network, data, scripts, results.

## Important files:

Data/accessories_v2.rda  -- an R data file with the following objects:
  - binMap -- bins to categorize samples by cell type score abundance
  - ccs -- the cell-cell scaffold (ccs)
  - cellECDFS -- the cell score abundance distributions

Data/cell_type_expr/cell_type_expression_tables.rda
  - The source data that was used to create the cell-cell scaffold.
  
Data/cell_type_expr/new_Celltype_map_march8.tsv
  - The map of cell names in source data to xCell data.


## Table of contents:

1. Analysis - Data and notebooks used in the analysis of networks.

2. Code - the code to build the cell-cell scaffold, weight the networks, and perform permuation testing.

3. Data - objects needed for computing the edge weights.

4. Expr - columns cut from the TCGA Pan-Cancer batch corrected RNA-seq data.

5. Perm - example of the permutation outputs.

6. Prob - cell abundance probability distribution examples.

7. Results - example result files showing edge weights.

8. SQL - BigQuery SQLs used to compute statistics.

9. testEnv - a complete test example to run through.


