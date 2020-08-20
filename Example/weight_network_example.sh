
# running a small portion of the cytokine network
# on a few samples

# need R packages:
#install.packages('hash')
#install.packages('data.table')

Rscript Code/compute_weights_v2_example.R Data/SampleExpr/sample_1234.rda 1234
Rscript Code/compute_weights_v2_example.R Data/SampleExpr/sample_42.rda 42
Rscript Code/compute_weights_v2_example.R Data/SampleExpr/sample_5403.rda 5403
Rscript Code/compute_weights_v2_example.R Data/SampleExpr/sample_6701.rda 6701

# run by chmod +x weight_network_example.sh
# and ./Examples/weight_network_example.sh

# the results are now in Example/Results/

# and in R:
#
#> load('Example/Results/res_1234.rda')
#> ls()
#[1] "res1" "si"   "ti"  
#> res1
#$Barcode
#[1] "TCGA-BA-4076-01A-01R-1436-07"
#
#$Vals
#[1] 0.170252127 0.006095121 0.068846329 0.221304214 0.104754050
#
#> si
#[1] "Data/SampleExpr/sample_1234.rda"
#> ti
#[1] "1234"


# ## checking the results ## #
#### against the BigQuery results database ####
#WITH
#  edgeinfo AS (
#  SELECT
#    EdgeID
#  FROM
#    `isb-cgc-02-0001.Cytokine_Network_Work.Cytokine_Network_Scaffold`
#  WHERE
#    Lcell = 'DC'
#    AND Ligand = 'CALM1'
#    AND Receptor = 'GRM4'
#    AND RCell = 'CD8+ T-cells' )
#    
#SELECT
#  *
#FROM
#  `isb-cgc-02-0001.Cytokine_Network_Work.weighted_cytokine_network`
#WHERE
#   EdgeID = (select EdgeID from edgeinfo) AND
#   AliquotBarcode = "TCGA-BA-4076-01A-01R-1436-07"
 
#### RESULTS FROM BIGQUERY ###
#Query complete (1.3 sec elapsed, 493.3 GB processed)
#Row	AliquotBarcode	SampleID	EdgeID	EdgeWt	
#1	TCGA-BA-4076-01A-01R-1436-07  1234   295806   0.1702521272

# confirmed!
   