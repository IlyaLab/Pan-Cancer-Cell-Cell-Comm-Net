#!/usr/bin/python
import sys
import pandas
import numpy as np

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


def featureGen(x, featureValue):
	# if x is the feature value, then true, else false
	if x == featureValue:
		return('1')
	else:
		return('0')


def idxLookup(d,b):
	if b in d:
		return(d[b])
	else:
		return(-1)


def main():
	# aliquotPath is path to ebpp barcode file
	# tablePath is path to the data table containing feature
	# featureName is found in the first row of the data table
	# barcodeName is name of column containing barcodes in table
	# barcodeMode: is it patientBarcode, sampleBarcode, aliquotBarcode
	# Rmode: past the first row, there's an extra row name to discard
	if (len(sys.argv) < 6):
		print('aliquotPath, tablePath, featureName, featureValue, barcodeName, barcodeMode, Rmode') 
		return(0)

	print('running')
	aliquotPath = sys.argv[1]
	tablePath   = sys.argv[2]
	featureName = sys.argv[3]
	featureValue= sys.argv[4]
	barcodeName = sys.argv[5]
	barcodeMode = sys.argv[6]
	Rmode       = sys.argv[7]
	foutPath    = sys.argv[8]

	ebpp = open(aliquotPath,'r').read().strip().split('\n')
	ebppFormat = formatBarcodes(ebpp, barcodeMode)
	ebppDict = dict([(xi,i) for i,xi in enumerate(ebppFormat)]) # barcode will give idx
	
	# read the clinical table
	print('reading' + tablePath + '\n')
	clin = pandas.read_csv(tablePath)

	# get the coded phenotype #
	pheno = [ featureGen(x, featureValue) for x in clin[featureName] ]

	# get the clin table barcodes #
	bcode = clin[barcodeName]
	bidx  = [idxLookup(ebppDict,bi) for bi in bcode]

	print('number of true cases:')
	print(np.sum([pi == '1' for pi in pheno]))
	print('number of samples:')
	print(len(pheno))
	print(len(bidx))

	fout = open(foutPath,'w')
	fout.write('\t'.join([str(b) for b in bidx ]) +'\n')
	fout.write('\t'.join([str(p) for p in pheno]) +'\n')
	fout.close()
	return(0)



if __name__ == "__main__":
    main()
