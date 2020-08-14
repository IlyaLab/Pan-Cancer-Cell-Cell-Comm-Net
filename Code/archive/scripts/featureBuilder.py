#!/usr/bin/python
import sys

# given ebpp aliquot barcodes in a file
# and a data table with column names that include
# barcodes and  some feature of interest.

# goal: produce (1) feature for ML and (2) the map of those 
# features to the ebpp columns (starting at 2? first column 
# has gene IDs.

def formatBarcodes(ebpp, mode):

	if mode == 'aliquots':
		return(ebpp)
	if mode == 'samples':
		return([x[0:16] for x in ebpp])
	if mode == 'patients':
		return([x[0:12] for x in ebpp])
	else:
		return([])


def main():
	# aliquotPath is path to ebpp barcode file
	# tablePath is path to the data table containing feature
	# featureName is found in the first row of the data table
	# barcodeName is name of column containing barcodes in table
	# barcodeMode: is it patientBarcode, sampleBarcode, aliquotBarcode
	# Rmode: past the first row, there's an extra row name to discard
	if (len(sys.argv) < 6):
		print('aliquotPath, tablePath, featureName, barcodeName, barcodeMode, Rmode') 
		return(0)

	print('running')
	aliquotPath = sys.argv[1]
	tablePath   = sys.argv[2]
	featureName = sys.argv[3]
	barcodeName = sys.argv[4]
	barcodeMode = sys.argv[5]
	Rmode       = sys.argv[6]

	ebpp = open(aliquotPath,'r').read().strip().split('\n')
	ebppFormat = formatBarcodes(ebpp, barcodeMode)

	


	return(0)



if __name__ == "__main__":
    main()
