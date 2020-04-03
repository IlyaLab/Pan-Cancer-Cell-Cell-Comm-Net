#aliquotPath, tablePath, featureName, featureValue, barcodeName, barcodeMode, Rmode
#bcr_patient_barcode,type,PFI.1
python3 featureBuilder_1vsall.py PanCancer/EBPP_Aliquots.csv PanCancer/PanCancer_Clinical_Extra_Endpoints.csv type PAAD bcr_patient_barcode patients False paad_vs_all_pheno.tsv
